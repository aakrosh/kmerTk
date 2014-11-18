#include "fastq.h"

// open the fastq file corresponding to the name of the file
FastqSequence* OpenFastqSequence(const char* const file,
                                 const Bool is_illumina_encoded,
                                 const Bool do_trim_reads) {
    FastqSequence* fqs = CkalloczOrDie(sizeof(struct FastqSequence_st));

    if (SameString(file, "-") == TRUE) {
        fqs->fp = stdin;
    } else {
        fqs->fp = CkopenOrDie(file, "r");
    }

    fqs->name_length = 1;
    fqs->name = CkalloczOrDie(fqs->name_length + 1);

    fqs->bases_length = 1;
    fqs->bases = CkalloczOrDie(fqs->bases_length + 1);

    fqs->quals_length = 1;
    fqs->quals = CkalloczOrDie(fqs->quals_length + 1);

    fqs->is_illumina_encoded = is_illumina_encoded;
    fqs->do_trim = do_trim_reads;

    return fqs;
}

//  free the memory allocated to the fqSequence structure and make it NULL
static void FreeSequence(FastqSequence** pfqSequence) {
    FastqSequence* sp = *pfqSequence;
    if (sp == NULL) return;
    Ckfree(sp->name);
    Ckfree(sp->bases);
    Ckfree(sp->quals);
    if (sp->fp != stdin) fclose(sp->fp);
    Ckfree(sp);
}

// read the next fastq sequence. Return FALSE when you get to the end of the
// file
static Bool ReadNextSequence(FastqSequence* sp) {
    size_t idx;

    // is this the end of the file
    ssize_t num_read;
    if ((num_read = Getline(&sp->name, &sp->name_length, sp->fp)) == -1) {
        FreeSequence(&sp);
        return FALSE;
    }
    if (num_read == 0) {
        FreeSequence(&sp);
        return FALSE;
    }
    ForceAssert(sp->name[0] == '@');
    idx = num_read - 1;
    while (sp->name[idx] != '\n') {
        sp->name[idx] = 0;
        idx--;
    }
    ForceAssert(sp->name[idx] == '\n');
    sp->name[idx] = 0;

    if ((num_read = Getline(&sp->bases, &sp->bases_length, sp->fp)) == -1) {
        PrintMessageThenDie("expected the bases for the read %s\n", sp->name);
    }
    idx = num_read - 1;
    while (sp->bases[idx] != '\n') {
        sp->bases[idx] = 0;
        idx--;
    }
    ForceAssert(sp->bases[idx] == '\n');
    sp->bases[idx] = 0;
    size_t sequence_length = idx;

    size_t  length_tmp = 1;
    char* qual_tmp = CkalloczOrDie(length_tmp + 1);
    if (Getline(&qual_tmp, &length_tmp, sp->fp) == -1) {
        PrintMessageThenDie("expected + and quals for the read %s\n", sp->name);
    }
    ForceAssert(qual_tmp[0] == '+');
    Ckfree(qual_tmp);

    if ((num_read = Getline(&sp->quals, &sp->quals_length, sp->fp)) == -1) {
        PrintMessageThenDie("expected + and quals for the read %s\n", sp->name);
    }
    idx = num_read - 1;
    while (idx >= sequence_length) {
        sp->quals[idx] = 0;  
        idx--;
    }

    // change the encoding if required
    if (sp->is_illumina_encoded == TRUE) {
        for (idx = 0; idx < sp->quals_length - 1; idx++) {
            sp->quals[idx] -= 31;
        }
    }

    // do we need to trim the reads
    idx = sp->quals_length - 1;
    sp->bases[idx] = 0;
    sp->quals[idx] = 0;
    if (sp->do_trim == TRUE) {
        while (idx > 0) {
            --idx;
            if ((sp->quals[idx] - 33) <= 2) {
                sp->bases[idx] = 0;
                sp->quals[idx] = 0;
            }
        }
    }

    sp->slen = idx;
    return TRUE; 
}


// open the fastq file and return the first FastqSequence
FastqSequence* ReadFastqSequence(const char* const file,
                                 const Bool is_illumina_encoded,
                                 const Bool do_trim) {
    FastqSequence* sp = OpenFastqSequence(file, is_illumina_encoded, do_trim);

    if (ReadNextSequence(sp) == FALSE) {
        return NULL;
    }

    return sp;
}

// get the next fastq FastqSequence to the FastqSequence sp
FastqSequence* GetNextSequence(FastqSequence* const sp) {
    Bool is_another_sequence = ReadNextSequence(sp);

    return is_another_sequence == TRUE ? sp : NULL;
}

// print the fastq FastqSequence
void PrintFastqSequence(const FastqSequence* const sp) {
    printf("%s\n", sp->name);
    printf("%s\n", sp->bases);
    printf("+\n");
    printf("%s\n", sp->quals);
}

// free all the resources used by this fastq FastqSequence
void CloseFastqSequence(FastqSequence* sp) {
    FreeSequence(&sp);  
}
