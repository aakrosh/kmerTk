#include <iostream>

#ifdef __cplusplus 
extern "C" {
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#endif
#include "inttypes.h"

#include "utilities.h"
#include "clparsing.h"
#include "kmer.h"
#include "fastq.h"
#ifdef __cplusplus
}
#endif

#include "sparse_hash.h"

/*
 * The pipeline is as follows:
 * 1. Read and store all the non-error kmers seen in the homogametic sex. 
 * 2. Read and store all the non-error kmers seen in the heterogametic sex.
 * 3. Iterate through the reads from the heterogametic individual, and tag
 *    every read as coming from one of the three:
 *      a) autosomes shared by both sexes.
 *      b) sex chromosome shared by both sexes.
 *      c) sex chromosome unique to the heterogametic sex.
 * 
 *  Step 1,2 will be done externally using either jellyfish or kmerz_count. The
 *  code in this source file is for step 3.
 *
 * NUANCES
 * -------
 * a) If we use PE sequences, the tagging will be done per fragment, so that 
 *    this information can be used downstream in assemblies for example.
 * b) Instead of using the reads, some might want to use an assembly, and that
 *    will be supported as well.
 */

Bool debug_flag;

const char* program_name = "kmerz_tag";
const char* program_description = 
"find novel sequences in one set from another set";
const char* program_use = 
"kmerz_tag [options] kmer_length kmers1 kmers2 fq2_1.txt [fq2_2.txt]";
const char* program_version_major = "0";
const char* program_version_minor = "1";
const char* program_revision_date = "11252014";

int main(int argc UNUSED, char** argv)
{
    // start book keeping
    t0 = time(0);

    // parse the command line
    CommandLineArguments* cl_options = NewCommandLineArguments();

    ParseOptions(&cl_options, &argc, &argv);

    // does the user just want some help
    Bool print_help = GetOptionBoolValueOrDie(cl_options, "help");
    if (print_help == TRUE) {
        PrintSimpleUsageString(cl_options);
        return EXIT_SUCCESS;
    }

    // does the user know the correct syntax
    if ((argc != 5) && (argc != 6)) {
        PrintSimpleUsageString(cl_options);
        return EXIT_FAILURE;
    }

    uint kmer_length;
    if (sscanf(argv[1], "%u", &kmer_length) != 1) {
        #ifdef Large
        PrintMessageThenDie("Kmer length should be an odd integer < 64: %s",
        argv[1]);
        #else
        PrintMessageThenDie("Kmer length should be an odd integer < 32: %s",
        argv[1]);
        #endif
    }

    char* hom_kmers_name = argv[2];
    char* het_kmers_name = argv[3];
    char* het_fq1_name   = argv[4];
    char* het_fq2_name   = NULL;
    if (argc == 6) {
        het_fq2_name = argv[5];
    }

    FILE* fp = NULL;
    size_t n = 1;
    char* fptr = (char*)CkallocOrDie(n+1);
    char* kmer_buffer = (char*)CkallocOrDie(1024);
    Kmer word, antiword, stored;
    uint64_t num_kmers_added = 0;
    int count;

    // add the kmers from the homogametic sex to this sparse hash set
    SparseHashMap hom_kmers;
    fp = CkopenOrDie(hom_kmers_name, "r");
    num_kmers_added = 0;
    while (Getline(&fptr, &n, fp) != -1) {
        if (sscanf(fptr, "%s %d\n", kmer_buffer, &count) != 2) {
            PrintMessageThenDie("Unable to parse %s", fptr);
        }

        word = BuildIndex(kmer_buffer, kmer_length);
        word = GetNextKmer(word, kmer_buffer, kmer_length, 0);
        antiword = ReverseComplementKmer(word, kmer_length);
        stored = word < antiword ? word : antiword;
        
        if (CheckKmerInSparseHashMap(hom_kmers, stored) == FALSE) {
            hom_kmers[stored] = count;
            num_kmers_added += 1;
        }
    }
    fclose(fp);
    fprintf(stderr, 
    "Added %"PRIu64" kmers from the homogametic sex\n", num_kmers_added);
    
    // now iterate through the heterogametic sex kmers to find the kmers seen
    SparseHashMap het_kmers;
    fp = CkopenOrDie(het_kmers_name, "r");
    num_kmers_added = 0;
    while (Getline(&fptr, &n, fp) != -1) {
        if (sscanf(fptr, "%s %d\n", kmer_buffer, &count) != 2) {
            PrintMessageThenDie("Unable to parse %s", fptr);
        }

        word = BuildIndex(kmer_buffer, kmer_length);
        word = GetNextKmer(word, kmer_buffer, kmer_length, 0);
        antiword = ReverseComplementKmer(word, kmer_length);
        stored = word < antiword ? word : antiword;
        
        if (CheckKmerInSparseHashMap(het_kmers, stored) == FALSE) {
            het_kmers[stored] = count;
            num_kmers_added += 1;
        }
    }
    fclose(fp);
    Ckfree(fptr);
    fprintf(stderr, 
    "Added %"PRIu64" kmers from the heterogametic sex\n", num_kmers_added);
 
    // now iterate through the fastq sequence from the heterogametic sex and
    // find ones that could be from the sex chromosome unique to this sex.
    FastqSequence* sequence1 = NULL;  
    FastqSequence* sequence2 = NULL;  
    sequence1 = ReadFastqSequence(het_fq1_name, FALSE, FALSE);
    if (het_fq2_name != NULL) {
        sequence2 = ReadFastqSequence(het_fq2_name, FALSE, FALSE);
    }

    while (sequence1) {
        uint num_kmers  = 0; // kmers that are not errors
        uint num_novel  = 0; // kmers only seen in the heterogametic sex
        uint num_common = 0; // kmers that are seen in both sexes

        // total number of kmers I will likely encounter in this fragment. This
        // includes all the error kmers, and assumes that all kmers in this
        // fragment are unique
        uint total_kmers = sequence1->slen - kmer_length + 1;
        if (sequence2 != NULL) {
            total_kmers += (sequence2->slen - kmer_length + 1);
        }        
        uint hom_coverage[total_kmers];
        uint het_coverage[total_kmers];
        uint indx = 0;

        if (sequence1->slen >= kmer_length) {
            // let account for all the kmers in this sequence
            word = BuildIndex(sequence1->bases, kmer_length);

            for (uint i = 0; i <= (sequence1->slen - kmer_length); i++) {
                word = GetNextKmer(word, sequence1->bases, kmer_length, i);
                antiword = ReverseComplementKmer(word, kmer_length);
                stored = word < antiword ? word : antiword;

                Bool het = CheckKmerInSparseHashMap(het_kmers, stored);
                Bool hom = CheckKmerInSparseHashMap(hom_kmers, stored);

                if ((het == FALSE) && (hom == FALSE)) {
                    // definitely an error kmer; dont count it
                    hom_coverage[indx] = 0;
                    het_coverage[indx] = 0;
                } else if ((het == TRUE) && (hom == FALSE)) {
                    // only seen in the heterogamtic sex; definitely interesting
                    num_novel += 1;
                    num_kmers += 1;
                    hom_coverage[indx] = 0;
                    het_coverage[indx] = het_kmers[stored];
                } else if ((het == FALSE) && (hom == TRUE)) {
                    // only seen in homogametic sex; probably a polymorphism?
                    num_kmers += 1;
                    hom_coverage[indx] = hom_kmers[stored];
                    het_coverage[indx] = 0;
                } else if ((het == TRUE) && (hom == TRUE)) {
                    // seen in both sexes
                    num_common += 1;
                    num_kmers += 1;
                    hom_coverage[indx] = hom_kmers[stored];
                    het_coverage[indx] = het_kmers[stored];
                } else {
                    PrintThenDie("not possible");
                }
                indx += 1;
            }
        }

        if ((sequence2 != NULL) && (sequence2->slen >= kmer_length)) {
            // let account for all the kmers in this sequence
            word = BuildIndex(sequence2->bases, kmer_length);

            for (uint i = 0; i <= (sequence2->slen - kmer_length); i++) {
                word = GetNextKmer(word, sequence2->bases, kmer_length, i);
                antiword = ReverseComplementKmer(word, kmer_length);
                stored = word < antiword ? word : antiword;

                Bool het = CheckKmerInSparseHashMap(het_kmers, stored);
                Bool hom = CheckKmerInSparseHashMap(hom_kmers, stored);

                if ((het == FALSE) && (hom == FALSE)) {
                    // definitely an error kmer; dont count it
                    hom_coverage[indx] = 0;
                    het_coverage[indx] = 0;
                } else if ((het == TRUE) && (hom == FALSE)) {
                    // only seen in the heterogamtic sex; definitely interesting
                    num_novel += 1;
                    num_kmers += 1;
                    hom_coverage[indx] = 0;
                    het_coverage[indx] = het_kmers[stored];
                } else if ((het == FALSE) && (hom == TRUE)) {
                    // only seen in homogametic sex; probably a polymorphism?
                    num_kmers += 1;
                    hom_coverage[indx] = hom_kmers[stored];
                    het_coverage[indx] = 0;
                } else if ((het == TRUE) && (hom == TRUE)) {
                    // seen in both sexes
                    num_common += 1;
                    num_kmers += 1;
                    hom_coverage[indx] = hom_kmers[stored];
                    het_coverage[indx] = het_kmers[stored];
                } else {
                    PrintThenDie("not possible");
                }
                indx += 1;
            }
        }

        // if the fraction of novel kmers is extremely high, then in most
        // probability this is a read that originates from the heterogametic
        // specific sex chromosome. If the fraction of common kmers is extremely
        // high, then it could either be from the autosomes, or the sex
        // chromosome that is shared by the two sexes.
        if (num_kmers > 0) {
            float frac_novel = num_novel * 1.0 / num_kmers;
            float frac_common = num_common * 1.0 / num_kmers;

            if (frac_novel == 1.0) {
                fprintf(stderr, "%s : %u : %u : %u -> chrY\n", sequence1->name, num_kmers, num_common, num_novel);
            } else if (frac_common == 1.0) {
                for (uint i = 0; i < indx; i++) {
                    fprintf(stderr, "%u\t", hom_coverage[i]);
                }               
                fprintf(stderr, "\n");
                for (uint i = 0; i < indx; i++) {
                    fprintf(stderr, "%u\t", het_coverage[i]);
                }               
                fprintf(stderr, "\n");
                fprintf(stderr, "%s : %u : %u : %u -> common\n", sequence1->name, num_kmers, num_common, num_novel);
            }
        }
        fprintf(stderr, "%s : %u : %u : %u\n", sequence1->name, num_kmers, num_common, num_novel);

        sequence1 = GetNextSequence(sequence1);
        if (sequence2) {
            sequence2 = GetNextSequence(sequence2);
        }
    }

    Ckfree(kmer_buffer);
    return EXIT_SUCCESS;
}
