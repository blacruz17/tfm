#!/bin/bash

##### metaWRAP #####
# RAW_READS should contain uncompressed fastq files
# for all the paired-end reads. Filename must be
# RUNACCESSION_1.fastq or RUNACCESSION_2.fastq

conda activate metaWRAP

mkdir READ_QC
for F in RAW_READS/*_1.fastq; do 
	R=${F%_*}_2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	metawrap read_qc -1 $F -2 $R -t 6 -o READ_QC/$SAMPLE
done > metawrapOUT.txt 2>&1 # obtains log file


##### concatenating runs #####
cd READ_QC
mkdir concatenated_runs
./concatenator.sh to-concatenate.txt concatenated-runs/

# to-concatenate.txt is a text file 
# each line contains runs that should be concatenated together
# separated by semicolons
# note: this script works for samples with 2 runs per sample
#		if there are extra runs, add for_n and rev_n
#		and do not forget to add them in your text file too

# concatenated-runs/ is the directory where 
# output files will be created

##### move clean reads without host ######
mkdir CLEAN_READS

for i in READ_QC/*; do 
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done

##### same for concatenated runs: #####
# move directory to CLEAN_READS
for f in READ_QC/concatenated-runs/; do 
	mv $f CLEAN_READS/
done

# move files to CLEAN_READS and
# remove the empty directory:
cd CLEAN_READS/concatenated-runs
mv * ..
cd ..
rm -r concatenated-runs

# fix filenames:
# this is necessary before running the metaphlan.pl script
cd ..
conc_files=$(ls CLEAN_READS/ | grep -E '_._..fastq')
for file in $conc_files; do
    new_name=$(sed 's/\_//' <(echo $file))
    mv CLEAN_READS/$file CLEAN_READS/$new_name
done

##### running MetaPhlAn #####
cd CLEAN_READS
conda activate metaphlan
# the only argument for metaphlan-script.pl is a list of
# forward reads ending with _1 (reverse shoud end with _2)
# metaphlanOUT.txt will be our log file
./metaphlan-script.pl *_1.fastq > metaphlanOUT.txt 2>&1

# after running metaphlan, merge tables:
## NOTE: run this command in the directory 
## containing the metaphlan output files
merge_metaphlan_tables.py profiled_*.txt > merged_abundance_table.txt