BINARIES = kmerz_count kmerz_select_enriched

CC  = gcc
CC1 = g++

CFLAGS  = -W -Wformat -Wimplicit -Wreturn-type -Wall -Werror \
          -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
          -Wcast-align -Winline -Wextra 
CPFLAGS = -W -Wformat -Wreturn-type -Wall -w \
          -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
          -Wcast-align -Winline -Wextra 
VERSION_FLAGS= \
	-DVERSION_MAJOR="\"${VERSION_MAJOR}"\" \
    -DVERSION_MINOR="\"${VERSION_MINOR}"\" \
    -DREVISION_DATE="\"${REVISION_DATE}"\" \

CFLAGS += -g -ggdb
CFLAGS += ${VERSION_FLAGS} 
CFLAGS += -O3

CPFLAGS += -g -ggdb
CPFLAGS += ${VERSION_FLAGS}
CPFLAGS += -O3 

SPINC = /afs/bx.psu.edu/home/ratan/include

# compile the version allowing upto 32mers
all: $(BINARIES)

# compile the version allowing upto 64mers
large: CFLAGS  += -DLarge 
large: CPFLAGS += -DLarge 
large: $(BINARIES)

.PHONY: clean

clean:
	-rm $(BINARIES) *.o
	-rm -rf *.dSYM

kmerz_count: utilities.h utilities.c sllist.h sllist.c clparsing.h clparsing.c kmer.h kmer.c murmur_hash.h murmur_hash.c fasta.h fasta.c fastq.h fastq.c kmerz_count.c
	$(CC)  $(CFLAGS) -c utilities.c
	$(CC)  $(CFLAGS) -c sllist.c
	$(CC)  $(CFLAGS) -c clparsing.c
	$(CC)  $(CFLAGS) -c fasta.c
	$(CC)  $(CFLAGS) -c fastq.c
	$(CC)  $(CFLAGS) -c kmer.c
	$(CC)  $(CFLAGS) -c murmur_hash.c
	$(CC)  $(CFLAGS) -c bloom_filter.c
	$(CC1) $(CPFLAGS) -I$(SPINC) -fopenmp -o kmerz_count utilities.o sllist.o clparsing.o kmer.o murmur_hash.o fasta.o fastq.o bloom_filter.o kmerz_count.c -lz -lm

kmerz_select_enriched: utilities.h utilities.c sllist.h sllist.c clparsing.h clparsing.c kmer.h kmer.c murmur_hash.h murmur_hash.c fasta.h fasta.c fastq.h fastq.c kmerz_select_enriched.c
	$(CC)  $(CFLAGS) -c utilities.c
	$(CC)  $(CFLAGS) -c sllist.c
	$(CC)  $(CFLAGS) -c clparsing.c
	$(CC)  $(CFLAGS) -c fasta.c
	$(CC)  $(CFLAGS) -c fastq.c
	$(CC)  $(CFLAGS) -c kmer.c
	$(CC)  $(CFLAGS) -c murmur_hash.c
	$(CC)  $(CFLAGS) -c bloom_filter.c
	$(CC1) $(CPFLAGS) -I$(SPINC) -fopenmp -o kmerz_select_enriched utilities.o sllist.o clparsing.o kmer.o murmur_hash.o fasta.o fastq.o bloom_filter.o kmerz_select_enriched.c -lz -lm


