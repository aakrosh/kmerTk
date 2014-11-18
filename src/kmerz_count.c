// This program is used to count the kmers (upto 64mers) in a fasta/q dataset

extern "C" {
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include "stdlib.h"
#include "inttypes.h"
#include "omp.h"

#include "clparsing.h"
#include "kmer.h"
#include "fasta.h"
#include "fastq.h"
#include "bloom_filter.h"
}

#include "sparse_hash.h"

#define chunk 1000000

const char* program_version_major = "0";
const char* program_version_minor = "1";
const char* program_revision_date = "11172014";
const char* program_name          = "kmerz_count";
const char* program_description   = 
    "count kmers in the input fastq/fasta dataset";
const char* program_use           =
"kmerz_count [options] kmer_length max_kmers_exp genome_size file.fq/file.fa";

Bool debug_flag;

static void ProcessChunk(char** const chunkbases,
                         const int numseqs,
                         const uint kmer_length,
                         SparseHashMap& kmers)
{
    int indx, tid;
    Kmer word, antiword, stored;

    std::vector<Kmer> array;
    std::vector<Kmer>::iterator it;
#pragma omp for
    for (indx = 0; indx < numseqs; indx += 1) {
        // let account for all the kmers in this sequence
        word = BuildIndex(chunkbases[indx], kmer_length);
        uint num_kmers = strlen(chunkbases[indx]) - kmer_length;

        for (uint i = 0; i <= num_kmers; i++) {
            word = GetNextKmer(word, chunkbases[indx],kmer_length,i);
            antiword = ReverseComplementKmer(word, kmer_length);
            stored = word < antiword ? word : antiword;
            
            array.push_back(stored);
        }

        Ckfree(chunkbases[indx]);
    }

#pragma omp critical
    for (it = array.begin(); it != array.end(); ++it) {
        if (CheckKmerInSparseHashMap(kmers, (*it)) == TRUE) {
            if (kmers[(*it)] <= (umaxof(Kcount) - 1)) {
                kmers[(*it)] += 1;
            }
        } else {
            kmers[(*it)] = 1;
        }
    }
}


static void ProcessChunkFirstRound(char** const chunkbases,
                                   const int numseqs,
                                   const uint kmer_length,
                                   BloomFilter* const bf,
                                   SparseHashMap& kmers)
{
    int indx, tid;
    Kmer word, antiword, stored;

    std::vector<Kmer> array;
    std::vector<Kmer>::iterator it;
#pragma omp for
    for (indx = 0; indx < numseqs; indx += 1) {
        // let account for all the kmers in this sequence
        word = BuildIndex(chunkbases[indx], kmer_length);
        uint num_kmers = strlen(chunkbases[indx]) - kmer_length;

        for (uint i = 0; i <= num_kmers; i++) {
            word = GetNextKmer(word, chunkbases[indx],kmer_length,i);
            antiword = ReverseComplementKmer(word, kmer_length);
            stored = word < antiword ? word : antiword;

            if (CheckKmerInSparseHashMap(kmers, stored) == FALSE) {
                if (CheckKmerInBloomFilter(bf, stored) == TRUE) {
                    array.push_back(stored);
                } else {
                    if (!AddKmerToBloomFilter(bf, stored)) {
                        array.push_back(stored);
                    }
                }
            } else {
                array.push_back(stored);
            }
        }

        Ckfree(chunkbases[indx]);
    }

#pragma omp critical
    for (it = array.begin(); it != array.end(); ++it) {
        kmers[(*it)] = 0;
    }
}

static void ProcessChunkSecondRound(char** const chunkbases,
                                    const int numseqs,
                                    const uint kmer_length,
                                    SparseHashMap& kmers)
{
    int indx, tid;
    Kmer word, antiword, stored;

    std::vector<Kmer> array;
    std::vector<Kmer>::iterator it;
#pragma omp for
    for (indx = 0; indx < numseqs; indx += 1) {
        // let account for all the kmers in this sequence
        word = BuildIndex(chunkbases[indx], kmer_length);
        uint num_kmers = strlen(chunkbases[indx]) - kmer_length;

        for (uint i = 0; i <= num_kmers; i++) {
            word = GetNextKmer(word, chunkbases[indx],kmer_length,i);
            antiword = ReverseComplementKmer(word, kmer_length);
            stored = word < antiword ? word : antiword;

            if (CheckKmerInSparseHashMap(kmers, stored) == TRUE) {
                array.push_back(stored);
            }
        }

        Ckfree(chunkbases[indx]);
    }

#pragma omp critical
    for (it = array.begin(); it != array.end(); ++it) {
        if (kmers[(*it)] <= (umaxof(Kcount) - 1)) {
            kmers[(*it)] += 1;
        }
    }
}

// count the kmers in the input file. Print the counts once you are done.
static void CountKmers(char* const* const file_names,
                       const int num_files,
                       const uint kmer_length,
                       const uint64_t num_expected_kmers,
                       const uint64_t expected_genome_size,
                       const int ploidy,
                       const Bool is_illumina_encoded,
                       const Bool do_trim,
                       const int progress_chunk,
                       const Bool is_fastq,
                       const int num_threads,
                       const Bool count_all) {
    omp_set_num_threads(num_threads);

    // this is the bloom filter to tag all the kmers that have been seen at
    // least once.
    BloomFilter* bf = NULL;
    if (count_all == FALSE) {
        bf = NewBloomFilter(0.001, num_expected_kmers, 0);
    }

    // this is the hash table I am going to use to store the counts
    SparseHashMap kmers;
    kmers.rehash(ploidy*expected_genome_size);

    // gather a chunk of sequences so that we can use multiple threads to deal
    // with them.
    char** chunkbases = (char**)CkallocOrDie(chunk * sizeof(char*));
    char* kmer_buffer = (char*)CkalloczOrDie(kmer_length + 1);

    void* sequence;
    uint64_t num_kmers_added = 0;
    uint64_t num_sequence_processed = 0;

    // go through all the files once to find the kmers
    for (int findx = 0; findx < num_files; findx++) {
        char* file_name = file_names[findx];

        if (is_fastq) {
            sequence = ReadFastqSequence(file_name,
                                         is_illumina_encoded,
                                         do_trim);
        } else {
            sequence = ReadFastaSequence(file_name);
        }
    
        uint bases_indx = 0;
    
        while (sequence) {
            uint slen;
            char* name;
            char* bases;
            if (is_fastq) {
                slen = ((FastqSequence*)sequence)->slen;
                name = ((FastqSequence*)sequence)->name;
                bases= ((FastqSequence*)sequence)->bases;
            } else {
                slen = ((FastaSequence*)sequence)->slen;
                name = ((FastaSequence*)sequence)->header;
                bases= ((FastaSequence*)sequence)->bases;
            }
    
            if (slen >= kmer_length) {
                if (debug_flag == TRUE) {
                    PrintDebugMessage("1. Processing %s", name);
                }
                // print progress
                num_sequence_processed += 1;
                if ((num_sequence_processed - 1) % progress_chunk == 0) {
                    double frac = 0.0;
                    if (bf) frac = bf->num_set_bits * 1.0 / bf->num_bits;
                    PrintDebugMessage(
                    "1. Processing read %"PRIu64": %s, kmers: %zu, bf :%.5f",
                    num_sequence_processed, is_fastq ? name + 1 : name, 
                    kmers.size(), frac);
                }
    
                // a load factor greater than 0.7-0.8 is a sign that the user did
                // not select the expected number of kmers judiciously. Lets warn
                // the user, as increasing the size of the hashtable can be very
                // slow.
                if (kmers.load_factor() > 0.8) {
                    PrintWarning("Current load factor: %2.6f",  
                    kmers.load_factor());
                    PrintWarning(
                    "Maybe try increasing expected genome size from %"PRIu64, 
                    expected_genome_size);
                }
    
                chunkbases[bases_indx++] = CopyString(bases);
                if ((bases_indx % chunk) == 0) {
#pragma omp parallel
                    {
                    if (count_all == TRUE) {
                        ProcessChunk(chunkbases,chunk,kmer_length,kmers);   
                    } else {
                        ProcessChunkFirstRound(chunkbases,chunk,kmer_length,bf,kmers);
                    }
                    }
                    bases_indx = 0;
                }
            }
    
            if (is_fastq) {
                sequence = GetNextSequence((FastqSequence*)sequence);
            } else {
                sequence = GetNextFastaSequence((FastaSequence*)sequence);
            }
        }
    
#pragma omp parallel
        {
        if (count_all == TRUE) {
            ProcessChunk(chunkbases, bases_indx, kmer_length, kmers);
        } else {
            ProcessChunkFirstRound(chunkbases, bases_indx, kmer_length, bf, kmers);
        }
        }

        PrintDebugMessage("1. Done with all the sequences in %s", file_name);
    }

    PrintDebugMessage("1. Counted %zu different kmers", kmers.size());
    ReportMemoryUsage();

    if (is_fastq) CloseFastqSequence((FastqSequence*)sequence);
    else CloseFastaSequence((FastaSequence*)sequence);
    
    if (count_all == TRUE) {
        goto printkmers;
    } else {
        FreeBloomFilter(&bf);
    }

    num_sequence_processed = 0;

    // lets count a second time and update the counts
    for (int findx = 0; findx < num_files; findx++) {
        char* file_name = file_names[findx];
        PrintDebugMessage("Reopening  %s to extract correct counts", file_name);

        if (is_fastq) {
            sequence = ReadFastqSequence(file_name,
                                         is_illumina_encoded,
                                         do_trim);
        } else {
            sequence = ReadFastaSequence(file_name);
        }
    
        uint bases_indx = 0;
        while (sequence) {
            uint slen;
            char* name;
            char* bases;
            if (is_fastq) {
                slen = ((FastqSequence*)sequence)->slen;
                name = ((FastqSequence*)sequence)->name;
                bases= ((FastqSequence*)sequence)->bases;
            } else {
                slen = ((FastaSequence*)sequence)->slen;
                name = ((FastaSequence*)sequence)->header;
                bases= ((FastaSequence*)sequence)->bases;
            }
    
            if (slen >= kmer_length) {
                if (debug_flag == TRUE) {
                    PrintDebugMessage("2. Processing %s", name);
                }
                // print progress
                num_sequence_processed += 1;
                if ((num_sequence_processed - 1) % progress_chunk == 0) {
                    double frac = bf->num_set_bits * 1.0 / bf->num_bits;
                    PrintDebugMessage(
                    "2. Processing read %"PRIu64": %s",
                    num_sequence_processed, is_fastq ? name + 1 : name);
                }
    
                chunkbases[bases_indx++] = CopyString(bases);
                if ((bases_indx % chunk) == 0) {
#pragma omp parallel
                    {
                    ProcessChunkSecondRound(chunkbases,chunk,kmer_length,kmers);
                    }
                    bases_indx = 0;
                }
            }
    
            if (is_fastq) {
                sequence = GetNextSequence((FastqSequence*)sequence);
            } else {
                sequence = GetNextFastaSequence((FastaSequence*)sequence);
            }
        }
    
#pragma omp parallel
        {
        ProcessChunkSecondRound(chunkbases, bases_indx, kmer_length, kmers);
        }
    }
    
    // now lets print the kmer counts
printkmers:
    SparseHashMap::const_iterator iter;
    for (iter = kmers.begin(); iter != kmers.end(); iter++) {
        const Kmer tmp = (*iter).first;
        ConvertKmerToString(tmp, kmer_length, &kmer_buffer);
        Kcount val = (*iter).second;
        if (sizeof(Kcount) == 1) {  
            if (count_all == TRUE) {
                printf("%s %"PRIu8"\n", kmer_buffer, val);
            } else {
                if (val > 1) printf("%s %"PRIu8"\n", kmer_buffer, val);
            }
        } else if (sizeof(Kcount) == 2) {
            if (count_all == TRUE) {
                printf("%s %"PRIu16"\n", kmer_buffer, val);
            } else {
                if (val > 1) printf("%s %"PRIu16"\n", kmer_buffer, val);
            }
        } else {
            PrintThenDie("This use-case has not been handled yet.");
        }
    }

    Ckfree(chunkbases);
    Ckfree(kmer_buffer);
}

int main(int argc, char** argv) {
    // start clock book-keeping 
    t0 = time(0);

    // parse the command line
    CommandLineArguments* cl_options = NewCommandLineArguments();

    // these are the valid options for the various commands
    AddOption(&cl_options, "threads", "1", TRUE, TRUE,
    "number of threads", NULL);
    AddOption(&cl_options, "format", "fastq", TRUE, TRUE,
    "format of input file (fastq/fasta)", NULL);
    AddOption(&cl_options, "qformat", "sanger", TRUE, TRUE,
    "encoding of quality (illumina/sanger).", NULL);
    AddOption(&cl_options, "trim", "FALSE", FALSE, TRUE,
    "should I trim the low quality 3' ends of the reads.", NULL);
    AddOption(&cl_options, "progress", "100000", TRUE, TRUE,
    "print progress every so many sequences", NULL);
    AddOption(&cl_options, "ploidy", "2", TRUE, TRUE,
    "ploidy for the species", NULL);
    AddOption(&cl_options, "exact", "FALSE", FALSE, TRUE,
    "count all kmers, including singletons",    NULL);
    
    ParseOptions(&cl_options, &argc, &argv);

    // does the user just want some help
    Bool print_help = GetOptionBoolValueOrDie(cl_options, "help");
    if (print_help == TRUE) {
        PrintSimpleUsageString(cl_options);
        return EXIT_SUCCESS;
    }

    // does the user know the correct syntax
    if (argc < 5) {
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
    if (kmer_length % 2 == 0) {
        PrintWarning("Kmer length should be an odd integer, using %u",
        --kmer_length);
    }

    uint64_t num_expected_kmers;
    if (sscanf(argv[2], "%"PRIu64, &num_expected_kmers) != 1) {
        PrintMessageThenDie("Number of expected kmers should be an integer: %s",
        argv[2]);
    }

    uint64_t expected_genome_size;
    if (sscanf(argv[3], "%"PRIu64, &expected_genome_size) != 1) {
        PrintMessageThenDie("Expected genome size should be an integer: %s",
        argv[3]);
    }

    // do I need additional debug info
    debug_flag = GetOptionBoolValueOrDie(cl_options, "debug");

    // count the kmers in these sequences
    Bool is_illumina_encoded = SameString(GetOptionStringValue(cl_options,
                                          "qformat"),"illumina")
                             ? TRUE : FALSE;
    Bool do_trim = GetOptionBoolValueOrDie(cl_options, "trim");
    int progress_chunk = GetOptionIntValueOrDie(cl_options, "progress");
    Bool is_fastq = TRUE;
    if (SameString(GetOptionStringValue(cl_options, "format"), "fasta")) {
        is_fastq = FALSE;
    }
    int num_threads = GetOptionIntValueOrDie(cl_options, "threads");
    Bool count_all = GetOptionBoolValueOrDie(cl_options, "exact");
    int ploidy = GetOptionIntValueOrDie(cl_options, "ploidy");

    //char* file_name = argv[4];
    CountKmers(argv + 4,
               argc - 4,
               kmer_length,
               num_expected_kmers,
               expected_genome_size,         
               ploidy,
               is_illumina_encoded,
               do_trim,
               progress_chunk,
               is_fastq,
               num_threads,
               count_all);

    FreeParseOptions(&cl_options, &argv);
    return EXIT_SUCCESS;
}
