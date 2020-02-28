#!/usr/bin/perl
use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
open PHEN, "All_samples_Exome_QC.phen";
my $queue="shortq"; # shortq: 240 hour, longq 480 hour
my $ppn=1;
while(<PHEN>){
chomp;
if(/FID/){
my $i;
my @phen=split(/\t/);
foreach my $phen(@phen){
next if $phen=~/FID/;
next if $phen=~/IID/;
$i++;
print $i;

print "$phen\n";
my $job_file_name = $phen . ".job";
my $status_file = $phen.".status";
my $curr_dir = $dir;
    open(OUT, ">$job_file_name") || die("Error in opening file $job_file_name.\n");
    print OUT "#!/bin/csh\n";
    print OUT "#PBS -N $phen\n";
    print OUT "#PBS -q $queue\n";  # glean is free, pdafm
    print OUT "#PBS -l nodes=1:ppn=$ppn\n";
    #print OUT "#PBS -o ".$phen.".log\n";
    #print OUT "#PBS -e ".$phen.".err\n";
    print OUT "#PBS -V\n";
    #print OUT "#PBS -M Guo.Shicheng\@marshfieldresearch.org \n";
    print OUT "#PBS -m abe\n";
    print OUT "cd $curr_dir\n";
    #print OUT "plink --bfile All_samples_Exome_QC --ci 0.95 --genotypic --remove excludeSample.txt --maf 0.01 --pheno All_samples_Exome_QC.phen --mpheno $i --covar All_samples_Exome_QC.cov --linear --out $phen[$i+1]\n";    
    #print OUT "plink --bfile All_samples_Exome_QC --linear mperm=50000 --ci 0.95 --seed 6377474 --remove excludeSample.txt --maf 0.01 --pheno All_samples_Exome_QC.phen --mpheno $i --covar All_samples_Exome_QC.cov --out $phen[$i+1]\n";    
    print OUT "plink --bfile All_samples_Exome_QC --linear mperm=50000  --covar All_samples_Exome_QC.cov --allow-no-sex --ci 0.95 --seed 6377474 --remove excludeSample.txt --maf 0.01 --pheno All_samples_Exome_QC.phen --mpheno $i --out $phen[$i+1]\n";       
    #print OUT "plink --bfile All_samples_Exome_QC --linear mperm=50000 --allow-no-sex --ci 0.95 --seed 6377474 --remove excludeSample.txt --maf 0.01 --pheno All_samples_Exome_QC.phen --mpheno $i --out $phen[$i+1].nonage\n";       
    close(OUT);
    system("qsub $job_file_name");
    }
  }
}



