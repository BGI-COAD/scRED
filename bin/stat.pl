#!usr/bin/perl -w
use strict;
die "perl $0 Tumor.anno.editing.txt Tumor.TypeDist.stat.txt Tumor.GenomeDist.stat.txt\n" unless (@ARGV==3);
my ($in,$out1,$out2) = @ARGV;
if ($in =~ /gz$/) {open IN, "<:gzip",$in or die $!;}
else {open IN, $in or die $!;}
open IN, $in or die $!;
open OUT1, ">$out1" or die $!;
open OUT2, ">$out2" or die $!;
my $numline = 0;
my $numlineall = 0;
my %hashv;
my %hashd;
<IN>;
while (<IN>){
    chomp;
    my $region;
    my @array = split /\s+/,$_;
    if ($array[-1] =~ /_/){$region = (split /_/,$array[-1])[1];}
    else {$region = $array[-1];}
    $numlineall += 1;
    if (exists $hashd{$array[-2]}){
        $hashd{$array[-2]} += 1;
    }
    else {
        $hashd{$array[-2]} = 1;
    }
    if (exists $hashd{$region}){
        $hashd{$region} += 1;
    }
    else {
        $hashd{$region} = 1;
    }
    if ($array[7] ne "\."){
        $numline += 1;
        if (exists $hashv{$array[7]}){
            $hashv{$array[7]} += 1;
        }
        else {
            $hashv{$array[7]} = 1;
        }
    }
}
close IN;
foreach (sort keys %hashv){
    my $p = $hashv{$_}/$numline;
    print OUT1 $_."\t".$hashv{$_}."\t".$p."\n"
}
close OUT1;
foreach (sort keys %hashd){
    my $p = $hashd{$_}/$numlineall;
    print OUT2 $_."\t".$hashd{$_}."\t".$p."\n"
}
close OUT2;

