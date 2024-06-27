####################### Example library #############
## Data analysis pipeline part 2 for SAM using UMI clustering. 
## Copyright © 2024 The Trustees of Columbia University in the City of New York. All Rights Reserved.


#!/bin/bash
#!/usr/bin/perl

## Replace the {PATH to the fastq.gz files} to the actual path
cd {PATH to the fastq.gz files} 
pwd
source /miniconda3/etc/profile.d/conda.sh

conda activate longread_umi
conda deactivate

parallel


#time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; gzip *.fastq; cd ..; pwd; done

####### end of skip  ############################

###### need to use the sed under longread_umi environment ######
 
conda activate longread_umi
conda info

 ############ retrieve the ori sequences for UMIs. 90% similarity #######
 ########## MS method #########
 ###### the final cluster count in the title of centroid fasta file 
 

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.fasta -fastaout test.umi12.f.u.fasta -sizeout -minuniquesize 1 -relabel umi -strand both -tabbedout test.umi12.f.u.tabbed; mkdir cluster_90; usearch -cluster_fast test.umi12.f.u.fasta -id 0.90 -centroids test.umi12.f.u.c.fasta -uc test.umi12.f.u.c.txt -sizein -sizeout  -strand both -minsize 1 -clusters cluster_90/test.umi12.f.u.c.clusters; grep ">" test.umi12.f.u.c.fasta | sed 's/^>//g' > test.umi12.f.u.c.title; cut -d ";" -f2 test.umi12.f.u.c.title | sort | uniq -c > test.umi12.f.u.c.title.count; cd ..; pwd;  done

######## result #########


######### end of results #############

###############################################################################
 ############ retrieve the ori sequences for UMIs. 94% similarity 2/36 #######
 ##### assuming the deunique is done at the 90% UMI #######

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -fastx_uniques test.umi12.f.fasta -fastaout test.umi12.f.u.fasta -sizeout -minuniquesize 1 -relabel umi -strand both -tabbedout test.umi12.f.u.tabbed; cd ..; pwd;  done


################
########### end of results ##########


time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_94; usearch -cluster_fast test.umi12.f.u.fasta -id 0.94 -centroids test.umi12.f.u.c.94.fasta -uc test.umi12.f.u.c.94.txt -sizein -sizeout -strand both -minsize 1 -sort size -clusters cluster_94/test.umi12.f.u.c.94.clusters; grep ">" test.umi12.f.u.c.94.fasta | sed 's/^>//g' > test.umi12.f.u.c.94.title; cut -d ";" -f2 test.umi12.f.u.c.94.title | sort | uniq -c > test.umi12.f.u.c.94.title.count; cd ..; pwd; done

################# result  

############ end of result  
#################################################################
############   make hist (histogram) file  ######################

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd
	
	awk -F ";|=" '{print $3;}' test.umi12.f.u.c.94.title > test.umi12.f.u.c.94.title.hist.txt
	wc -l test.umi12.f.u.c.94.title.hist.txt
 
	cd ..
	pwd

done


################### First analyze clusters with 8 or more reads ##########################

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; mkdir cluster_94_8up; usearch -cluster_fast test.umi12.f.u.fasta -id 0.94 -centroids test.umi12.f.u.c.94.8up.fasta -uc test.umi12.f.u.c.94.8up.txt -sizein -sizeout  -strand both -minsize 8 -sort size -clusters cluster_94_8up/test.umi12.f.u.c.94.8up.clusters; grep ">" test.umi12.f.u.c.94.8up.fasta | sed 's/^>//g' > test.umi12.f.u.c.94.8up.title; cut -d ";" -f1 test.umi12.f.u.c.94.8up.title > test.umi12.f.u.c.94.8up.title.ID; cut -d ";" -f2 test.umi12.f.u.c.94.8up.title | sort | uniq -c > test.umi12.f.u.c.94.8up.title.count; cd ..; pwd; done

################### result


#########################

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.u.c.94.8up.title.ID test.umi12.f.u.c.94.8up.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.u.c.94.8up.clusters; wc -l test.umi12.f.u.c.94.8up.title.ID; wc -l test.umi12.f.u.c.94.8up.clusters; cd ..; pwd; done


conda activate longread_umi

########### Example genomic location of mutation: 17: 24720177
################# extract sam file for each cluster, 8up  ######################
#########not the longread_umi samtools 
conda deactivate

echo "8up;"

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	sed 's/^/test.umi12.f.u.c.94.8up.clusters/g' test.umi12.f.u.c.94.8up.clusters > test.umi12.f.u.c.94.8up.clusters.rename
	 
	cd cluster_94_8up;
	pwd;
 
	time parallel -j 12 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.u.c.94.8up.clusters.rename;
 
	time parallel -j 12 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.cut.g.a.GRCh37.sam > {}.sam; samtools view -bT /Users/sw3203/Documents/Fun/Ref/hg19/Bowtie2Index/genome.fa {}.sam | samtools sort - {}; samtools index {}.bam; rm {}.sam; rm {}.seq;" < ../test.umi12.f.u.c.94.8up.clusters.rename;
 
	cd ../..
	pwd

done

###############  Calling variants in each cluster with DP>0.5 ####################
for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd
	
	time parallel -j 12 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f /Users/sw3203/Documents/Fun/Ref/hg19/Bowtie2Index/genome.fa {}.bam > {}.vcf ' < ../test.umi12.f.u.c.94.8up.clusters.rename;
		

	time parallel -j 12 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>4" {}.vcf > {}.50.4.vcf; bgzip -c {}.50.4.vcf > {}.50.4.vcf.gz; tabix {}.50.4.vcf.gz; grep -v "#" {}.50.4.vcf > {}.50.vcf ' < ../test.umi12.f.u.c.94.8up.clusters.rename;
	
	 
	 cd ../..
done

	

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd
	
	rm 50.vcf.count
	rm call.50.vcf
	touch 50.vcf.count
	touch call.50.vcf
	
	for f in $(cat ../test.umi12.f.u.c.94.8up.clusters.rename); do wc -l  $f.50.vcf >> 50.vcf.count; cat $f.50.vcf >> call.50.vcf;  done
	awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	cd ../..
	pwd
	
done

########### Example genomic location of mutation: 17: 24720177 

time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; cd cluster_94_8up; pwd; grep "256972" call.50.vcf | wc -l;  grep "256972" call.50.vcf | grep -v "DP=1;" | wc -l; cd ../..; pwd; done	

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd	
	awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	for f in $(cat 50.vcf.count.1); do echo $f; cat $f; done > 50.vcf.count.1.input
 ## check at 24720177, list all the supporting clusters
 ## Replace it with the actual genomic location of interest 
	grep -B1 "24720177" 50.vcf.count.1.input
	cut -f2 call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	
	tar -zcvf cluster_94.tar.gz cluster_94
	rm -r cluster_94

	cd ..
	pwd
	
done

###################################################################################
###################### analyze clusters with 4-7 reads ############################

conda activate longread_umi
conda info


########## use the cluster files in cluster_94_8up ###########


time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; awk -F ";|=" ' $3>3 && $3<8 {print $1}' test.umi12.f.u.c.94.title > test.umi12.f.u.c.94.4-7.title.ID; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.u.c.94.4-7.title.ID test.umi12.f.u.c.94.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.u.c.94.4-7.clusters; wc -l test.umi12.f.u.c.94.4-7.title.ID; wc -l test.umi12.f.u.c.94.4-7.clusters; sed 's/^/test.umi12.f.u.c.94.8up.clusters/g' test.umi12.f.u.c.94.4-7.clusters > test.umi12.f.u.c.94.4-7.clusters.rename; cd ..; pwd; done

########## results ###################


########## end of results ######################


################# sam for each cluster ##############
######### Don't use the samtools in the longread_umi   
conda deactivate


time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 	 
	cd cluster_94_8up;
	pwd;
 
	time parallel -j 12 ' grep "^>" {} | cut -d ";" -f1 | sed "s/^>//g" > {}.ID ' < ../test.umi12.f.u.c.94.4-7.clusters.rename;
 
	time parallel -j 12 " awk 'FNR==NR{a[\$0];next}(\$2 in a)' {}.ID ../test.umi12.f.u.tabbed > {}.seq; cut -f1 {}.seq > {}.seq.ID; awk 'FNR==NR{a[\$0];next}(\$1 in a)' {}.seq.ID ../${f%.fastq.gz}.cut.g.a.GRCh37.sam > {}.sam; samtools view -bT /Users/sw3203/Documents/Fun/Ref/hg19/Bowtie2Index/genome.fa {}.sam | samtools sort - {}; samtools index {}.bam; rm {}.sam; rm {}.seq;" < ../test.umi12.f.u.c.94.4-7.clusters.rename;
 
	cd ../..
	pwd

done



conda activate longread_umi

for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd
	
	time parallel -j 12 ' bcftools mpileup -Ov -d 10000 -L 10000 -a "FORMAT/AD,FORMAT/DP" -f /Users/sw3203/Documents/Fun/Ref/hg19/Bowtie2Index/genome.fa {}.bam > {}.vcf ' < ../test.umi12.f.u.c.94.4-7.clusters.rename;
		
	time parallel -j 12 'bcftools view -i "FORMAT/AD[0:1]/FORMAT/DP>0.5 && INFO/DP>4" {}.vcf > {}.50.4.vcf; bgzip -c {}.50.4.vcf > {}.50.4.vcf.gz; tabix {}.50.4.vcf.gz; grep -v "#" {}.50.4.vcf > {}.50.vcf ' < ../test.umi12.f.u.c.94.4-7.clusters.rename;
		 
	cd ../..
done
	

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd
	
	rm 4-7.50.vcf.count
	rm call.4-7.50.vcf
	touch 4-7.50.vcf.count
	touch call.4-7.50.vcf
	
	for f in $(cat ../test.umi12.f.u.c.94.4-7.clusters.rename); do wc -l  $f.50.vcf >> 4-7.50.vcf.count; cat $f.50.vcf >> call.4-7.50.vcf;  done
	awk '$1>0 {print $2;}' 4-7.50.vcf.count > 4-7.50.vcf.count.1
	cd ../..
	pwd
	
done


time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd	
#	awk '$1>0 {print $2;}' 4-7.50.vcf.count > 4-7.50.vcf.count.1
	for f in $(cat 4-7.50.vcf.count.1); do echo $f; cat $f; done > 4-7.50.vcf.count.1.input
	#grep -B1 "24720177" 50.vcf.count.1.input
	cut -f2 call.4-7.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done
   

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd	
	grep -B1 "24720177" 50.vcf.count.1.input
	grep -B1 "24720177" 4-7.50.vcf.count.1.input
	#cut -f2 call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;

	cat call.50.vcf call.4-7.50.vcf > call.4up.50.vcf
	echo "4up:"
	cut -f2 call.4up.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
#	pwd
	
done

########################################################################## 
############## cutoff at 94% avg-2 size ################

######### assuming the 8up is done ###################################
######### change the -minsize to the average coverage ###########

echo " 94% avg-2 size "
conda activate longread_umi


time for f in $(cat input.file); do echo $f; cd ${f%.fastq.gz}; pwd; usearch -cluster_fast test.umi12.f.u.fasta -id 0.94 -centroids test.umi12.f.u.c.94.Ave-2.fasta -uc test.umi12.f.u.c.94.Ave-2.txt -sizein -sizeout  -strand both -minsize 5 -sort size ; grep ">" test.umi12.f.u.c.94.Ave-2.fasta | sed 's/^>//g' > test.umi12.f.u.c.94.Ave-2.title; cut -d ";" -f1 test.umi12.f.u.c.94.Ave-2.title > test.umi12.f.u.c.94.Ave-2.title.ID; cut -d ";" -f2 test.umi12.f.u.c.94.Ave-2.title | sort | uniq -c > test.umi12.f.u.c.94.Ave-2.title.count; awk -F "\t|;" 'FNR==NR{a[$0];next}($9 in a)' test.umi12.f.u.c.94.Ave-2.title.ID test.umi12.f.u.c.94.Ave-2.txt | cut -f2 | sort -k1n | uniq > test.umi12.f.u.c.94.Ave-2.clusters; wc -l test.umi12.f.u.c.94.Ave-2.title.ID; wc -l test.umi12.f.u.c.94.Ave-2.clusters; sed 's/^/test.umi12.f.u.c.94.8up.clusters/g' test.umi12.f.u.c.94.Ave-2.clusters > test.umi12.f.u.c.94.Ave-2.clusters.rename; cd ..; pwd; done


   
time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	pwd 
	sed 's/^/test.umi12.f.u.c.94.8up.clusters/g' test.umi12.f.u.c.94.Ave-2.clusters > test.umi12.f.u.c.94.Ave-2.clusters.rename;
	cd ..
done
	
conda deactivate

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd
	
	rm Ave-2.50.vcf.count
	rm Ave-2.call.50.vcf
	touch Ave-2.50.vcf.count
	touch Ave-2.call.50.vcf
	
	for f in $(cat ../test.umi12.f.u.c.94.Ave-2.clusters.rename); do wc -l  $f.50.vcf >> Ave-2.50.vcf.count; cat $f.50.vcf >> Ave-2.call.50.vcf; done
	awk '$1>0 {print $2;}' Ave-2.50.vcf.count > Ave-2.50.vcf.count.1
	cd ../..
	pwd
	
done

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd	
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	for f in $(cat Ave-2.50.vcf.count.1); do echo $f; cat $f; done > Ave-2.50.vcf.count.1.input
	grep -B1 "24720177" 50.vcf.count.1.input
	cut -f2 Ave-2.call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	cd cluster_94_8up;
	pwd	
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	grep -B2 "24720177" Ave-2.50.vcf.count.1.input
	cut -f2 Ave-2.call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done


time for f in $(cat input.file)
do 
	echo $f; 
	cd ${f%.fastq.gz}; 
	 wc -l test.umi12.f.u.c.94.8up.title.ID; wc -l test.umi12.f.u.c.94.8up.clusters;
	cd cluster_94_8up;
	pwd	
	#awk '$1>0 {print $2;}' 50.vcf.count > 50.vcf.count.1
	#for f in $(cat 50.vcf.count.1); do echo $f; cat $f; done > 50.vcf.count.1.input
	grep -B1 "24720177" 50.vcf.count.1.input
	cut -f2 call.50.vcf | sort | uniq -c | sort -k1nr 
	cd ../..
	pwd
	
done

exit
