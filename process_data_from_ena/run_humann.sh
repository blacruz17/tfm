#!/bin/bash

##### concatenate files #####
# this script concatenates forward and reverse files
# must be run in the script containing the clean reads
# only argument is the list of forward files
cd CLEAN_READS/
perl cat_files.pl *1.fastq

##### run HUMAnN #####
## we can provide the directory containing the metaphlan
## taxonomic profiles to save some time:
## in this case profiles are in OUTPUT/profiles
## this loop removes temporary output while maintaining
## individual log files from HUMAnN for each sample
## a log for the whole process is also created
conda activate metaphlan
cd ..
for F in CLEAN_READS/*.fastq; do 
	BASE=${F##*/}
	TAX=OUTPUT/profiles/profiled_${BASE%.*}.txt
	humann -i $F -o HUMANN --verbose --taxonomic-profile $TAX
	rm HUMANN/${BASE%.*}_humann_temp/*bowtie2*
    rm HUMANN/${BASE%.*}_humann_temp/*aligned*
    rm HUMANN/${BASE%.*}_humann_temp/*phlan*
done > humannOUT.txt 2>&1


##### post-HUMAnN processing #####
# gets all files in the same directory:
for file in $(find -type f | grep tsv); do cp $file ../../humann_output_files/; done
# joins tables:
humann_join_tables -i hmp_subset -o hmp_subset_genefamilies.tsv --file_name genefamilies
humann_join_tables --input hmp_subset --output hmp_subset_pathcoverage.tsv --file_name pathcoverage
humann_join_tables --input hmp_subset --output hmp_subset_pathabundance.tsv --file_name pathabundance

## HUMAnN to STAMP:
humann_split_stratified_table -i all_pathway_abundances-relab.tsv -o split_folder
humann_split_stratified_table -i all_pathway_coverages.tsv -o split_folder_coverage

sed 's/pathway/MetaCyc_pathway/' split_folder/all_pathway_abundances-relab_unstratified.tsv > pathabundance_unstratified_head-fix.tsv
sed 's/pathway/MetaCyc_pathway/' split_folder_coverage/all_pathway_coverages_unstratified.tsv > pathcoverage_unstratified_head-fix.tsv

# fix headers:
sed -i 's/.conc_Abundance//g' pathabundance_unstratified_head-fix.tsv
sed -i 's/_Abundance//g' pathabundance_unstratified_head-fix.tsv

# change NAs to 0 and some other fixes:
sed -i 's/NA/0/g' pathabundance_unstratified_head-fix.tsv
sed -i 's/.conc_Coverage//g' pathcoverage_unstratified_head-fix.tsv
sed -i 's/_Coverage//g' pathcoverage_unstratified_head-fix.tsv
sed -i 's/NA/0/g' pathcoverage_unstratified_head-fix.tsv