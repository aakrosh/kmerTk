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
 * Read and store all the kmers seen in the homogametic sex. Then iterate
 * through the list of kmers from the heterogametic sex, and save the ones that
 * are unique to it. Traverse the reads or contigs of the heterogametic sex to
 * calculate whether what fraction of the kmers in that read/contigs are from
 * that novel set. 
 * 
 * NUANCES
 * a) We will use a sparseset to store the kmers.
 * b) We will use PE sequences, so the fraction will be calculated per read, but 
 *    if either of the reads qualify as described above, then I will print out 
 *    both the reads in separate files. This makes it convinient to use 
 *    in downstream assemblies
 * c) Instead of using the reads, some might want to use an assembly, and that 
 *    will be supported as well. In that case we will consider the number of
 *    novel kmers that we see in the assembly rather than the reads. 
 * 
 */

Bool debug_flag;

const char* program_name = "kmerz_select_enriched";
const char* program_description = 
"find novel sequences in one set from another set";
const char* program_use = 
"kmerz_select_enriched [options] kmer_length kmers1 kmers2 fq2.txt";
const char* program_version_major = "0";
const char* program_version_minor = "1";
const char* program_revision_date = "11062014";

int main(int argc UNUSED, char** argv)
{
    // start book keeping
    t0 = time(0);

    // parse the command line
    CommandLineArguments* cl_options = NewCommandLineArguments();

    AddOption(&cl_options, "fraction", "0.75", TRUE, TRUE,
    "fraction of kmers required to be unique", NULL);

    ParseOptions(&cl_options, &argc, &argv);

    // does the user just want some help
    Bool print_help = GetOptionBoolValueOrDie(cl_options, "help");
    if (print_help == TRUE) {
        PrintSimpleUsageString(cl_options);
        return EXIT_SUCCESS;
    }

    // does the user know the correct syntax
    if (argc != 5) {
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

    // what is the fraction of kmers that should be unique in a sequence for it
    // to be considered to be from the heterogametic sex.
    char* tmp = GetOptionStringValue(cl_options, "fraction");
    float fraction;
    if (sscanf(tmp, "%f", &fraction) != 1) {
        PrintMessageThenDie("Could not understand %s", tmp);
    }

    char* hom_kmers_name = argv[2];
    char* het_kmers_name = argv[3];
    char* het_fq_name    = argv[4];

    // add the kmers from the homogametic sex to this sparse hash set
    SparseHashSet hom_kmers;
    
    FILE* fp = CkopenOrDie(hom_kmers_name, "r");
    size_t n = 1;
    char* fptr = (char*)CkallocOrDie(n+1);

    char* kmer_buffer = (char*)CkallocOrDie(1024);
    Kmer word, antiword, stored;
    uint64_t num_kmers_added = 0;

    while (Getline(&fptr, &n, fp) != -1) {
        if (sscanf(fptr, "%s %*d\n", kmer_buffer) != 1) {
            PrintMessageThenDie("Unable to parse %s", fptr);
        }

        word = BuildIndex(kmer_buffer, kmer_length);

        word = GetNextKmer(word, kmer_buffer, kmer_length, 0);
        antiword = ReverseComplementKmer(word, kmer_length);
        stored = word < antiword ? word : antiword;
        
        if (CheckKmerInSparseHashSet(hom_kmers, stored) == FALSE) {
            hom_kmers.insert(stored);
            num_kmers_added += 1;
        }
    }
    fclose(fp);
    fprintf(stderr, 
    "Added %"PRIu64" kmers from the homogametic sex\n", num_kmers_added);
    
    // now iterate through the heterogametic sex kmers to find the kmers seen
    SparseHashSet het_kmers;

    fp = CkopenOrDie(het_kmers_name, "r");

    num_kmers_added = 0;
    while (Getline(&fptr, &n, fp) != -1) {
        if (sscanf(fptr, "%s %*d\n", kmer_buffer) != 1) {
            PrintMessageThenDie("Unable to parse %s", fptr);
        }

        word = BuildIndex(kmer_buffer, kmer_length);

        word = GetNextKmer(word, kmer_buffer, kmer_length, 0);
        antiword = ReverseComplementKmer(word, kmer_length);
        stored = word < antiword ? word : antiword;
        
        if (CheckKmerInSparseHashSet(het_kmers, stored) == FALSE) {
            het_kmers.insert(stored);
            num_kmers_added += 1;
        }
    }
    fclose(fp);
    Ckfree(fptr);
    fprintf(stderr, 
    "Added %"PRIu64" kmers from the heterogametic sex\n", num_kmers_added);
 
    // now iterate through the fastq sequence from the heterogametic sex and
    // find ones that could be from the sex chromosome unique to this sex.
    FastqSequence* sequence = ReadFastqSequence(het_fq_name, FALSE, FALSE);

    while (sequence) {
        if (sequence->slen >= kmer_length) {
            // let account for all the kmers in this sequence
            word = BuildIndex(sequence->bases, kmer_length);
            uint num_kmers = 0;  // kmers that are not errors
            uint num_novel = 0;  // kmers only seen in the heterogametic sex

            for (uint i = 0; i <= num_kmers; i++) {
                word = GetNextKmer(word, sequence->bases, kmer_length, i);
                antiword = ReverseComplementKmer(word, kmer_length);
                stored = word < antiword ? word : antiword;

                Bool het = CheckKmerInSparseHashSet(het_kmers, stored);
                Bool hom = CheckKmerInSparseHashSet(hom_kmers, stored);

                if ((het == FALSE) && (hom == FALSE)) {
                    // definitely an error kmer; dont count it
                } else if ((het == TRUE) && (hom == FALSE)) {
                    // only seen in the heterogamtic sex; definitely interesting
                    num_novel += 1;
                    num_kmers += 1;
                } else if ((het == FALSE) && (hom == TRUE)) {
                    // only seen in homogametic sex; probably a polymorphism?
                    
                } else if ((het == TRUE) && (hom == TRUE)) {
                    // seen in both sexes
                    num_kmers += 1;
                } else {
                    PrintThenDie("not possible");
                }
            }

            if (num_kmers > 0) {
                float frac = num_novel * 1.0 / num_kmers;
                if (frac >= fraction) {
                    //fprintf(stderr, "%d %d %f\n", num_novel, num_kmers, frac);
                    PrintFastqSequence(sequence);
                }
            }
        }

        sequence = GetNextSequence(sequence);
    }

    Ckfree(kmer_buffer);
    return EXIT_SUCCESS;
}