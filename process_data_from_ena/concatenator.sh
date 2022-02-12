#!/bin/bash
# textfile containing runs to be concatenated in each line
# runs in each line are separated by semicolons:
file=$1 
# directory where output files will be created:
conc_dir=$2

# concatenates files and removes originals:
while IFS= read -r line; do
	base_name=$(echo $line | cut -d ';' -f1)
	for_1=$(echo $line | cut -d ';' -f1)/final_pure_reads_1.fastq
	for_2=$(echo $line | cut -d ';' -f2)/final_pure_reads_1.fastq

	rev_1=$(echo $line | cut -d ';' -f1)/final_pure_reads_2.fastq
	rev_2=$(echo $line | cut -d ';' -f2)/final_pure_reads_2.fastq
	
	cat $for_1 $for_2 > $conc_dir/$base_name.conc_1.fastq
	cat $rev_1 $rev_2 > $conc_dir/$base_name.conc_2.fastq
	rm $for_1 $for_2 $rev_1 $rev_2
done < $file