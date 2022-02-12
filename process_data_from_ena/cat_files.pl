#!/usr/bin/perl
foreach $ar (@ARGV)
{
# gets forward (_1) file:
@field = split (/_/, $ar);
$name = $ar; # full filename
$name =~ s/_1.fastq//; # filename before the underscore

# gets the correct reverse (_2) file:
$r1 = $ar;
$r2 = $ar;
$r2 =~ s/1.fastq/2.fastq/;

# concatenates files, removes originals and 
# changes whitespaces to underscores:
system("cat $r1 $r2 > $name.fastq; rm $r1 $r2; sed -i 's/ /_/' $name.fastq");
}
