#!/usr/bin/perl
use strict;
use Cwd;

my @file=glob("*assoc.linear");
foreach my $file(@file){
        print "$file\n";
	open OUT, ">$file.pvalue";
        open F,$file || die "cannot open $file\n";
	while(<F>){
        chomp;
        next if /NA/;
        next if /AvgAge/;
        my @line=split/\s+/;
        print OUT "$_\n" if $line[-1] <0.00001;
        }
}
