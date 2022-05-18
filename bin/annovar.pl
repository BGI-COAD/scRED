#!usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
die "perl $0 in annovar refgene aluBed out.txt!" unless (@ARGV==5);
my ($in,$annovar,$refgene,$aluBed,$out) = @ARGV;

`perl $Bin/format.to.anno.pl $in $in.anno`;
`perl $annovar/annotate_variation.pl --splicing_threshold 5 -buildver mm10 $in.anno $annovar/mousedb_mm10`;

open RG, $refgene or die $!;
my %gene;
while (<RG>){
     chomp;
     my @c=split /\s+/,$_;
     $gene{$c[-4]}=$c[3];
}
close RG;

`perl $Bin/pos_filterByRegion.pl $aluBed $in >$in.alu.txt`;
my %hash;
open ALU,"$in.alu.txt" or die $!;
while(<ALU>){
    chomp;
    if($_=~/#/){next;}
    my @c=split /\s+/,$_;
    $hash{$c[0]}{$c[1]} = "Alu";
}
close ALU;

open OUT, ">$out" or die $!;
print OUT "Chromosome\tPosition\tReference\tAlteration\tDepth\tReference_support_reads\tAlteration_support_reads\tGenotype\tGene\tStrand\tAlu\tRegion\n";
open IN, "$in.anno.variant_function" or die $!;
my %comp = ("A"=>"T","T"=>"A","C"=>"G","G"=>"C");
while (<IN>){
    chomp;
    my $alu;
    my @c = split /\s+/,$_;
    if (exists $hash{$c[2]}{$c[3]}){
	$alu = "Alu";
    }
    else {$alu = "non-Alu";}
    if ($c[0] =~ /intergenic/ || $c[0] =~ /upstream/ ||  $c[0] =~ /downstream/){
        print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t"."."."\t".$c[1]."\t"."."."\t".$alu."\t".$c[0]."\n";
        next;
    }
    my $gn = $c[1];
    my @list;
    if ($gn =~ /,/ && $gn !~ /;/ && $gn !~ /:/){@list = split /,/,$gn;}
    if ($gn =~ /;/ && $gn !~ /,/ && $gn !~ /:/){@list = split /;/,$gn;}
    if ($gn =~ /;/ && $gn =~ /,/ && $gn !~ /:/){
        my @list1 = split /,/,$gn; 
        foreach (@list1){
            my @list2 = split /;/,$_;
            foreach my $a (@list2){
                push @list, $a;
            }
        }
    }
    if (@list == 0){
        $gn = (split /\(/,$gn)[0];
 	    if ($gene{$gn} eq "+"){
	        my $type = $c[5]."->".$c[6];
	        print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t".$type."\t".$c[1]."\t"."+"."\t".$alu."\t".$c[0]."\n";
     	}
    	if ($gene{$gn} eq "-"){
	        my $type = $comp{$c[5]}."->".$comp{$c[6]};
	        print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t".$type."\t".$c[1]."\t"."-"."\t".$alu."\t".$c[0]."\n";
    	}
    }
    else {
    	my $len = @list;
    	my $sl;
    	my $plus = "+" x $len;
    	my $minus = "-" x $len;
    	foreach (@list){
	        $gn = (split /\(/,$_)[0];
            if (!defined $gene{$gn}) {print $gn."\n";next;}
	        $sl .= $gene{$gn};
    	}
    	if ($sl eq $plus){
	        my $type = $c[5]."->".$c[6];
	        print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t".$type."\t".$c[1]."\t"."+"."\t".$alu."\t".$c[0]."\n";
    	}
    	elsif ($sl eq $minus){
	        my $type = $comp{$c[5]}."->".$comp{$c[6]};
	        print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t".$type."\t".$c[1]."\t"."-"."\t".$alu."\t".$c[0]."\n";
    	}
    	else {
	        if ($c[5] eq "T" && $c[6] eq "C"){
	        	my $type = $comp{$c[5]}."->".$comp{$c[6]};
	        	print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t".$type."\t".$c[1]."\t"."-"."\t".$alu."\t".$c[0]."\n";
            }
	        else {
                my $type = $c[5]."->".$c[6];
	        	print OUT $c[2]."\t".$c[3]."\t".$c[5]."\t".$c[6]."\t".$c[7]."\t".$c[8]."\t".$c[9]."\t".$type."\t".$c[1]."\t"."+"."\t".$alu."\t".$c[0]."\n";
            }
    	}
    }
}
close OUT;
