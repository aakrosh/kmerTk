/*! Routines for reading and writing fasta file */

#include "fasta.h"

/*all valid bases in the ASCII world should be 1*/
static const char bases[] = {
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,1,1,1,1,0,
0,1,1,0,0,1,0,1,1,0,
0,0,1,1,1,0,1,1,0,1,
0,0,0,0,0,0,0,1,1,1,
1,0,0,1,1,0,0,1,0,1,
1,0,0,0,1,1,1,0,1,1,
0,1,0,0,0,0,0,0,0
};


/*open the fasta file corresponding to the name of the file*/
FastaSequence* OpenFastaSequence(const char* const file)
{
	FastaSequence* sp = CkalloczOrDie(sizeof(struct FastaSequence_st));
	FILE* fp = CkopenOrDie(file,"r");

	sp->fp = fp;
	sp->rc = FALSE;
	return sp;
}

static void Append(char** parray,      /*the pointer to the array to be used*/
                   uint* const pmax,   /*how many bytes have been allocated*/
                   uint* const plen,   /*how many bytes have been used*/
                   const int ch)       /*the byte to be added*/
{
        uint max = *pmax;
        uint len = *plen;

        /*do we need more memory*/
        if(len >= max){
                *parray = CkreallocOrDie(*parray, max + CHUNK);
                *pmax = *pmax + CHUNK;
        }
        (*parray)[len] = ch;
        *plen = *plen + 1;
}

/*read the FastaSequence*/
static Bool _ReadFastaSequence(FastaSequence* const sp)
{
	/*is it the end of the file?*/
	if(feof(sp->fp) || ferror(sp->fp)){
		return FALSE;
	}
	
	int ch;

	/*skip all the whitespace characters*/
	do{
		ch = fgetc(sp->fp);
	}while(ch == ' ' || ch == '\t');

	/*check and add the fasta header*/
	assert(ch == '>');
	ch = fgetc(sp->fp);
	do{
		Append(&sp->header, &sp->hmax, &sp->hlen, ch);
		ch = fgetc(sp->fp);
	}while(ch != '\n');
	Append(&sp->header, &sp->hmax, &sp->hlen, 0);
	sp->hlen--;

	/*add the bases*/
	do{
		ch = fgetc(sp->fp);
		if(1 == bases[ch]){
			Append(&sp->bases, &sp->smax, &sp->slen, ch);
		}
	}while(ch != '>' && ch != EOF);
	Append(&sp->bases, &sp->smax, &sp->slen, 0);
	sp->slen--;

	if(ch != EOF && ungetc(ch, sp->fp) == EOF){
		PrintThenDie("error in pushing back the >");
	}
	
	return TRUE;
}

/*free the memory allocated to the FastaSequence structure and make it NULL*/
static void FreeSequence(FastaSequence** const pFastaSequence)
{
	if(*pFastaSequence){
		Ckfree((*pFastaSequence)->header);
		Ckfree((*pFastaSequence)->bases);
		if((*pFastaSequence)->fp != NULL && (*pFastaSequence)->fp != stdin){	
			fclose((*pFastaSequence)->fp);
		}
		Ckfree(*pFastaSequence);
		*pFastaSequence = NULL;
	}
}

/*open the fasta file and return the first FastaSequence*/
FastaSequence* ReadFastaSequence(const char* const file)
{
	/*open the fasta file*/
	FastaSequence* sp = OpenFastaSequence(file);

	/*read the FastaSequence*/
	if(_ReadFastaSequence(sp) == FALSE){
		return NULL;
	}

	return sp;
}

/*get the next fasta FastaSequence to the FastaSequence sp*/
FastaSequence* GetNextFastaSequence(FastaSequence* const sp)
{
	sp->slen = 0;
	sp->hlen = 0;
	if(_ReadFastaSequence(sp) == FALSE){
		return NULL;
	}

	return sp;
}

/*free the resources with this FastaSequence*/
void CloseFastaSequence(FastaSequence* sp)
{
	FreeSequence(&sp);
}

/*pretty print this FastaSequence*/
void Format(const char* const FastaSequence, const int len)
{
	int i,j;
	for(i = 0, j = 0; i < len; i++){
		if(j != 0  && (j % 60)==0){
			printf("\n");
		}
		if(FastaSequence[i] != '-'){
			printf("%c", FastaSequence[i]);
			j++;
		}
	}
	printf("\n");
}

/*print the fasta FastaSequence in a format with 60 bases in a single line, and 
 *with a fasta header*/
void PrintFastaSequence(const FastaSequence* const sp)
{
	printf(">%s\n", sp->header);
	Format(sp->bases, sp->slen);
}

/*clone the fasta FastaSequence*/
FastaSequence* CopySequence(const FastaSequence* const sp)
{
	FastaSequence* fsp = CkalloczOrDie(sizeof(struct FastaSequence_st));

	fsp->hmax = fsp->hlen = strlen((char*)sp->header);
	fsp->header = CopyString((char*)sp->header);

	fsp->smax = fsp->slen = strlen((char*)sp->bases);
	fsp->bases = CopyString((char*)sp->bases);
	return fsp;
}

/*reverse complement the FastaSequence in place*/
FastaSequence* ReverseComplementSequence(FastaSequence* const sp)
{
	sp->bases = ReverseComplementString(sp->bases, sp->slen);	
	sp->rc = sp->rc == TRUE ? FALSE : TRUE;
	return sp;
}

/*count the number of FastaSequences in the file*/
int CountSequences(const char* const file)
{
	FastaSequence* sp;
	int count = 0;

	if((sp = ReadFastaSequence(file)) == NULL){
		PrintMessageThenDie("error in reading the fasta reads from %s", file);
	}

	while(sp){
		count++;
		if(!GetNextFastaSequence(sp)){
			break;
		}
	}
	CloseFastaSequence(sp);

	return count;
}
