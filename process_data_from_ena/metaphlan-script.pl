#!/usr/bin/perl
foreach $ar (@ARGV)
{

@field = split (/_/, $ar);
# gets forward file:
$name = $ar; # name == full filename
$name =~ s/_1.fastq//; # name == filename up to underscore

# gets the correct reverse file:
$r1 = $ar; # r1 == full filename including _1
$r2 = $ar; # r2 == full filename including _1
$r2 =~ s/1.fastq/2.fastq/; # changes _1 to _2

# runs MetaPhlAn on both files:
system("metaphlan $r1,$r2 --bowtie2db ~/blanca/bowtiedb --bowtie2out $name.bowtie2.bz2 --nproc 6 --input_type fastq -o profiled_$name.txt");
# for PRJEB2054 we add --read_min_len 40,
# comment previous line and uncomment the next one:
# system("metaphlan $r1,$r2 --bowtie2db ~/blanca/bowtiedb --bowtie2out $name.bowtie2.bz2 --nproc 6 --input_type fastq -o profiled_$name.txt --read_min_len 40");
}