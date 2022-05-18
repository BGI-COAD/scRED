#!usr/bin/perl -w
use strict;
my ($in,$out) = @ARGV;
open IN, $in or die $!;
open OUT, ">$out" or die $!;
while (<IN>){
    chomp;
    my @array = split /\t/,$_;
    if ($array[10] eq "Alu"){
        print OUT $_."\n";
    }
    if ($array[10] eq "non-Alu" && $array[11] !~ /splicing/){
        print OUT $_."\n";
    }
}
close OUT;
close IN;
