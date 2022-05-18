#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
my $usage =<<USAGE;
Parameter

    --in           [STR]  a position file needed to filt
    --simpleRepeat [STR]  a simplerepeat locus file
    --reference    [STR]  a reference file
    --out          [STR]  the output file
    --endfre       [num]  the frequency of end for alt to filt [defaut 0.9]
    --endp         [num]  the p value of end fisher test to filt [defaut 0.05]
    --strandfre    [num]  the frequency of strand for alt to filt [defaut 0.9]
    --strandp      [num]  the p value of strandbias fisher test to filt [defaut 0.005]
    --homopolymer  [int]  the length of homopolymer [defaut 5]
    --help   give this information !

perl $0 --in postion.txt --simpleRepeat simpleRepeat.txt --reference reference.fa --out postion.filt.txt --endfre 0.9 --endp 0.05 --strandfre 0.9 --strandp 0.005 --homopolymer 5 

USAGE

my ($in,$simplerep,$ref,$out,$endfre,$endp,$strandfre,$strandp,$homopolymer,$help);
GetOptions (
    "in=s" =>\$in,
    "simpleRepeat=s"=>\$simplerep,
    "reference=s"=>\$ref,
    "out=s" =>\$out,
    "endfre=f" =>\$endfre,
    "endp=f" =>\$endp,
    "strandfre=f"=>\$strandfre,
    "strandp=f" =>\$strandp,
    "homopolymer=i" =>\$homopolymer,
    "help" =>\$help,
);
if (!$in || !$out || $help) {die $usage;}
if (!defined $endfre) {$endfre = 0.9;}
if (!defined $endp) {$endp = 0.05;}
if (!defined $strandfre) {$strandfre = 0.9;}
if (!defined $strandp) {$strandp = 0.005;}
if (!defined $homopolymer) {$homopolymer = 5;}
`perl $Bin/pos_filterByRegion.pl $simplerep $in > $out.rep`;
open REF, $ref or die $!;
$/ = ">";<REF>;
my %hashr;
while (<REF>){
    chomp;
    my @array = split /\n/,$_;
    my $chr = (split /\s+/,$array[0])[0];
    if ($chr !~ /chr/){$chr = "chr".$chr;}
    shift(@array);
    my $seq = join("",@array);
    $hashr{$chr} = $seq;
}
close REF;
$/ = "\n";
open SR, "$out.rep" or die $!;
my %hashs;
while (<SR>){
    chomp;
    my @array = split /\s+/,$_;
    my $pos = $array[0]."\t".$array[1];
    $hashs{$pos} = 1;
}
close SR;
`rm -f $out.rep`;
open IN,"gzip -dc $in |" or die $!;
open OUT,">$out" or die $!;
while(<IN>){
    chomp;
    if($_=~/^Chr/){
        print "$_\n";
    }else{
        my @c=split /\s+/,$_;
        my $chr = $c[0];
        if ($chr !~ /chr/){$chr = "chr".$c[0];}
        my $locus = $c[1];
        my $refbase = $c[3];
        my $pos = $c[0]."\t".$c[1];
        my $freq=$c[12]/$c[2];
        my $sb=$c[15]/($c[15]+$c[16]);
        my $refEndPor=0;
        my $altEndPor=0;
        my $poly = substr($hashr{$chr},$locus-1-$homopolymer+1,2*$homopolymer-1);
        if (exists $hashs{$pos}) {next;}
        if($c[12]<=2 || $c[12]<$c[19] || $freq<0.1 || $sb<0.1 || $sb>0.9){next;}
        if (($c[14]/($c[13]+$c[14]))>=$endfre){
            next;
        }
        if ($poly =~ /($refbase){$homopolymer,}/){next;}
        if($c[20]<$endp){
            $refEndPor=$c[6]/($c[5]+$c[6]) unless ($c[5]+$c[6]==0);
            $altEndPor=$c[14]/($c[13]+$c[14]) unless ($c[13]+$c[14]==0);
            if($altEndPor>$refEndPor){next;}
        }
        if ($c[16]/($c[15]+$c[16]) > $strandfre || $c[16]/($c[15]+$c[16]) < 1-$strandfre){
            next;
        }
        if($c[21]<$strandp){
            my $refstrandPor=0;
            $refstrandPor=abs(($c[8]/($c[7]+$c[8]))-0.5) unless ($c[7]+$c[8]==0);
            my $altstrandPor=0;
            $altstrandPor=abs(($c[16]/($c[15]+$c[16]))-0.5) unless ($c[15]+$c[16]==0);
            if($altstrandPor>$refstrandPor){next;}
        }
        print OUT "$_\n";
    }
}
close IN;
close OUT;
