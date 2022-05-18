#!usr/bin/perl -w
use strict;
my ($in,$out) = @ARGV;
open IN, $in or die $!;
open OUT, ">$out" or die $!;
while (<IN>){
    chomp;
    my @array = split /\s+/,$_;
	print OUT $array[0]."\t".$array[1]."\t".$array[1]."\t".$array[3]."\t".$array[11]."\t".$array[2]."\t".$array[4]."\t".$array[12]."\n";
}
close IN;
close OUT;
