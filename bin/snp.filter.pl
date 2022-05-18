#!/usr/bin/perl -w
use strict;
die "perl $0 input dbsnp output !\n" if @ARGV != 3;
my ($in,$dbsnp,$out) = @ARGV;

my %hash;
open IN, $in or die $!;
while (<IN>){
    chomp;
    my @array = split /\s+/,$_;
    my $chrpos = $array[0]."_".$array[1];
    $hash{$chrpos} = $_;
}
close IN;

if ($dbsnp =~ /\.gz$/){open DBSNP, "gzip -dc $dbsnp |" or die $!;}
else {open DBSNP, $dbsnp or die $!;}
while (<DBSNP>){
    chomp;
    if ($_ =~ /^#/) {next;}
    my @array = split /\s+/,$_;
    my $chrpos = "chr".$array[0]."_".$array[1];
    if (exists $hash{$chrpos}) {
        delete $hash{$chrpos};
    }
}
close DBSNP;

open OUT, ">$out" or die $!;
foreach my $each (keys %hash){
    print OUT $hash{$each}."\n";
}
close OUT;
