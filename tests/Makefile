SHELL := /bin/bash
C = ./../src
F = female_1.fq
T = 1

all: fetch truecount nonsinglecountfq nonsinglecountfa exactcountfq 

fetch: 
	wget http://www.bx.psu.edu/~ratan/kmerTk_data/male_1.fq
	wget http://www.bx.psu.edu/~ratan/kmerTk_data/male_2.fq
	wget http://www.bx.psu.edu/~ratan/kmerTk_data/female_1.fq
	wget http://www.bx.psu.edu/~ratan/kmerTk_data/female_2.fq

truecount: 
	@$C/simple_kmer_counter 31 $F > truecounts.txt

nonsinglecountfq: 
	@$C/kmerz_count --threads=$T 31 1150000 1130000 $F > nonsingletons.txt
	diff --brief <(cat truecounts.txt | awk '$$2 > 1 {print $$1,$$2}' | sort) <(cat nonsingletons.txt | sort)
	@rm nonsingletons.txt

nonsinglecountfa:
	@cat $F | awk '{if(NR%4==1) printf ">%s\n", substr($$1,2); else if(NR%4==2) print $$0}' > female.fa
	@$C/kmerz_count --threads=$T --format=fasta 31 1150000 1130000 female.fa > nonsingletons.txt
	diff --brief <(cat truecounts.txt | awk '$$2 > 1 {print $$1,$$2}' | sort) <(cat nonsingletons.txt | sort)
	@rm female.fa nonsingletons.txt

exactcountfq:
	@$C/kmerz_count --threads=$T --exact 31 1150000 1130000 $F > counts.txt
	diff --brief <(cat truecounts.txt | sort) <(cat counts.txt | sort)
	@rm counts.txt

findnovel:
	@$C/kmerz_count --threads=$T 31 2640000 140000 female_1.fq female_2.fq > female.kmers
	@$C/kmerz_count --threads=$T 31 3640000 180000 male_1.fq male_2.fq > male.kmers
	@$C/kmerz_select_enriched 31 female.kmers male.kmers <(cat male_1.fq male_2.fq) > male.only.reads
	@rm female.kmers male.kmers
	cat male.only.reads | awk 'NR % 4 == 1' | wc -l
	cat male.only.reads | awk 'NR % 4 == 1' | grep "sexchrom2" | wc -l
	@rm male.only.reads
