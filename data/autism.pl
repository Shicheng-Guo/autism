#!/usr/bin/perl
use strict;
open PHEN, "All_samples_Exome_QC.phen";
while(<PHEN>){
if(/FID/){
my $i;
my @phen=split(/\t/);
foreach my $phen(@phen){
next if $phen=~/FID/;
next if $phen=~/IID/;
$i++;
#system("plink --bfile All_samples_Exome_QC --ci 0.95 --genotypic --maf 0.01 --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno $i --covar All_samples_Exome_QC.cov --linear --out $phen[$i+1].ci");		
system("plink --bfile All_samples_Exome_QC --ci 0.95 --maf 0.01 --allow-no-sex --pheno All_samples_Exome_QC.phen --mpheno $i --covar All_samples_Exome_QC.cov --linear --out $phen[$i+1].ci");		
sleep(240)
}
}
}
#
