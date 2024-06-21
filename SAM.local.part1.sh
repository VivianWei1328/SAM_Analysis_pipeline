####################### SM-12 SPECC1L library #############
#!/bin/bash
#!/usr/bin/perl

cd /Users/sw3203/Documents/Research/Sequencing/illumina/UMI/SM-14_R2P
source /miniconda3/etc/profile.d/conda.sh

### need to install longread_umi conda environment as in README.md ####
conda activate longread_umi
conda deactivate

parallel

######################################################
#replace the reverse primer (R) with actual reverse primer (R), the reverse complement of forward primer (F-RC) with the actual F-RC 
# For example
#ATP6V0C.UMI-MR2T-200F	TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-NNNYRNNNYRNNNYRNNN-TCGTGGGCATGATCCTGATT
#ATP6V0C.UMI-MR1-200R	ACACTCTTTCCCTACACGACGCTCTTCCGATCT-NNNYRNNNYRNNNYRNNN-CAAGACCAACTACTGGGGGC
# R is CAAGACCAACTACTGGGGGC
# F-RC is AATCAGGATCATGCCCACGA
# 
# An example sequence from MiSeq is:
#>M04254:476:000000000-DK5K2:1:1101:14700:1774 1:N:0:1
#CTCCAAGGCACTCCGAGTCAAGACCAACTACTGGGGGCGGCGGCCCCGTGCGTATGTGTCAGGCTGTTCGTTCTGGAATGAGGAGGGGTGGTCTTTACATAATATTCTGTGGCTTGTGGGCTCGGAGAGGGTCTACTTTGTGGAGAGGATGAGGGCGACGATGAGACCGTAGAGGCCGAGCACCTCGGCGAAGATGAGAATCAGGATCATGCCCACGACTTTGAACCATAATATAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTCTGCACCTTCTCTTTTTCCTTCTTCTCTT
#
############ cut directly using cutadapt ########

############ replace -g with R primer, -a with RC of F primer ################

time for f in $(cat input.file); 
do echo $f; 
	cutadapt -g CAAGACCAACTACTGGGGGC -O 18 -e 0.05 -o ${f%.fastq.gz}.cut.g.fastq $f >${f%.fastq.gz}.cut.g.log; cutadapt -a AATCAGGATCATGCCCACGA -O 18 -e 0.05 -m 50 -o ${f%.fastq.gz}.cut.g.a.fastq ${f%.fastq.gz}.cut.g.fastq >${f%.fastq.gz}.cut.g.a.log; bowtie2 --threads 4 -x /Users/sw3203/Documents/Fun/Ref/hg19/Bowtie2Index/genome -U ${f%.fastq.gz}.cut.g.a.fastq -q --very-sensitive -S ${f%.fastq.gz}.cut.g.a.GRCh37.sam;
	samtools view -bS ${f%.fastq.gz}.cut.g.a.GRCh37.sam | samtools sort - ${f%.fastq.gz}.cut.g.a.GRCh37.sort; 
	samtools index ${f%.fastq.gz}.cut.g.a.GRCh37.sort.bam; 
	gzip $f;
	mkdir ${f%.fastq.gz};
	mv ${f%.fastq.gz}.* ${f%.fastq.gz};
done 

############ manually make u1 u2 fasta files ########
############ manually make u1 u2 fasta files ########

######### replace -a -g with R primer


######### replace -a -g with R primer
time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; cutadapt -a CAAGACCAACTACTGGGGGC -O 20 -e 0.05 -m 18 -M 18 --discard-untrimmed -o ${f%.fastq.gz}.u1.fastq ./$f >${f%.fastq.gz}.u1.log; cutadapt -g CAAGACCAACTACTGGGGGC -O 20 -e 0.05 --discard-untrimmed -o ${f%.fastq.gz}.u1.g.fastq $f >${f%.fastq.gz}.u1.g.log; cd ..; pwd; done

######### replace RC-F  (-g) with RC-F primer
time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; cutadapt -g AATCAGGATCATGCCCACGA -O 20 -e 0.05 --discard-untrimmed -o ${f%.fastq.gz}.u1.g.g.fastq ${f%.fastq.gz}.u1.g.fastq >${f%.fastq.gz}.u1.g.g.log; cutadapt -a AGATCGGAAGAGCACACGTCT -O 20 -e 0.05 -m 18 -M 18 --discard-untrimmed -o ${f%.fastq.gz}.u1.g.g.u2.fastq ${f%.fastq.gz}.u1.g.g.fastq >${f%.fastq.gz}.u1.g.g.u2.log; grep ^@ ${f%.fastq.gz}.u1.fastq | cut -d " " -f1 > ${f%.fastq.gz}.u1.id; cd ..; pwd; done

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; gunzip $f; awk 'FNR==NR{a[$0];next}($1 in a){print;getline;print;getline;print;getline;print;}' ${f%.fastq.gz}.u1.id ${f%.gz} > ${f%.fastq.gz}.u1.ori.fastq; grep ^@ ${f%.fastq.gz}.u1.g.g.u2.fastq | cut -d " " -f1 > ${f%.fastq.gz}.u1.g.g.u2.ID; awk 'FNR==NR{a[$0];next}($1 in a){print;getline;print;getline;print;getline;print;}' ${f%.fastq.gz}.u1.g.g.u2.ID ${f%.fastq.gz}.u1.fastq > ${f%.fastq.gz}.u1.u2pick.fastq; wc -l ${f%.fastq.gz}.u1.fastq; awk 'FNR==NR{a[$0];next}($1 in a){print;getline;print;getline;print;getline;print;}' ${f%.fastq.gz}.u1.id ${f%.fastq.gz}.u1.g.g.u2.fastq > ${f%.fastq.gz}.u2.u1.fastq; wc -l ${f%.fastq.gz}.u2.u1.fastq; grep "^@" ${f%.fastq.gz}.u2.u1.fastq | cut -d " " -f1 > ${f%.fastq.gz}.u2.u1.ID; awk 'FNR==NR{a[$0];next}($1 in a){print;getline;print;getline;print;getline;print;}' ${f%.fastq.gz}.u2.u1.ID ${f%.fastq.gz}.u1.fastq > ${f%.fastq.gz}.u1.u2.u1.fastq; wc -l ${f%.fastq.gz}.u1.u2.u1.fastq; gzip ${f%.gz}; cd ..; pwd; done


############### skip ##############################################
###### need to use the sed under longread_umi environment ######

echo "Run this part in terminal" 
cd /Users/sw3203/Documents/Research/Sequencing/illumina/UMI/SM-14_R2P
conda activate longread_umi
conda info
########## this step has syntex error in .sh file, need to work on ########

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; paste   <(sed -n '1~4s/^@/>/p;2~4p' ${f%.fastq.gz}.u1.u2.u1.fastq) <( sed -n '1~4s/^@/>/p;2~4p' ${f%.fastq.gz}.u2.u1.fastq ) |  cut -d " " -f1 | sed 's/\t//g' > ${f%.fastq.gz}.u12.fasta; wc -l ${f%.fastq.gz}.u12.fasta; PATTERN="[ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{6}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}[CT][AG][ATCG]{3}"; grep -B1 -E "$PATTERN" ${f%.fastq.gz}.u12.fasta | sed '/^--$/d' > ${f%.fastq.gz}.u12.f.fasta; wc -l ${f%.fastq.gz}.u12.f.fasta; cp ${f%.fastq.gz}.u12.f.fasta test.umi12.f.fasta; cd ..; pwd; done



time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; gzip *.fastq; cd ..; pwd; done

####### end of skip  ############################

















































