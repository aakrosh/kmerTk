/*! Routines for reading and writing fasta FastaSequences */

#ifndef FASTA_H_
#define FASTA_H_

#include "utilities.h"

/*bases in a single line of FastaSequences*/
#define FORMAT 60

#define CHUNK 200

typedef struct FastaSequence_st
{
	char* header;		/*the header of the FastaSequence*/
	uint hmax;			/*memory allocated for the header*/
	uint hlen;			/*length of the header*/

	char* bases;    	/*the base FastaSequence*/
	uint smax;			/*memory allocated for the FastaSequence*/
	uint slen;			/*length of the FastaSequence*/
	
	FILE* fp;			/*which file do we read*/
	Bool rc;			/*if set this FastaSequence is reverse complemented*/
}FastaSequence;

/*open the fasta file corresponding to the name of the file*/
FastaSequence* OpenFastaSequence(const char* const file);

/*open the fasta file and return the first FastaSequence*/
FastaSequence* ReadFastaSequence(const char* const file);

/*pretty print this FastaSequence*/
void Format(const char* const FastaSequence, const int len);

/*print the fasta FastaSequence in a format with 60 bases in a single line, and 
 *with a fasta header*/
void PrintFastaSequence(const FastaSequence* const sp);

/*get the next fasta FastaSequence to the FastaSequence sp*/
FastaSequence* GetNextFastaSequence(FastaSequence* const sp); 

/*free all the resources used by this fasta FastaSequence*/
void CloseFastaSequence(FastaSequence* const sp);

/*clone the fasta FastaSequence*/
FastaSequence* CopySequence(const FastaSequence* const sp);

/*reverse complement the FastaSequence in place*/
FastaSequence* ReverseComplementSequence(FastaSequence* const sp);

/*count the number of FastaSequences in the file*/
int CountSequences(const char* const file);

#endif
