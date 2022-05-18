#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
my $usage=<<USAGE;

Description:
    This program is used to identify and filter RNA editing sites based on tophat alignment file(bam). 

Parameters:
    --rnabam       [STR] the RNA alignment file which is used to calling RNA editing sites (required).
    --reference    [STR] the genome fasta file(hg19.fa) (required).
    --snpdb        [STR] a file including SNP database files (including directory and file), eg ./temp/dbsnp138, each row indicate a SNP database (required).
    --simpleRepeat [STR] genome simple repeat region annotation file, should be bed format (required).
    --alu          [STR] genome alu region annotation file, should be bed format (required).
    --dnabam       [STR] the DNA aligment file, it could be used to remove SNP [optional].
    --hg19_mRNA    [STR] a concat reference (combination hg19.fa and mRNA.fa), and variation reads would be mapped to this reference using bwa (required).
    --bwa          [STR] the directory of bwa (required).
    --samtools     [STR] the directory of samtools (required).
    --annovar      [STR] the directory of annovar (required).
    --output       [STR] the directory of output (required).
    --help         [STR] show this information!

Example:
    Perl $0 --rnabam rnabam --reference hg19.fa --dbsnp dbsnp.txt --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --hg19_mRNA hg_mrna.fa --annovar annovar --output test

USAGE

my ($rnabam,$reference,$snpdb,$simpleRepeat,$alu,$dnabam,$hg19_mRNA,$bwa,$samtools,$annovar,$output,$help);
GetOptions(
    "rnabam=s"=>\$rnabam,
    "reference=s"=>\$reference,
    "snpdb=s"=>\$snpdb,
    "simpleRepeat=s"=>\$simpleRepeat,
    "alu=s"=>\$alu,
    "dnabam=s"=>\$dnabam,
    "hg19_mRNA=s"=>\$hg19_mRNA,
    "bwa=s"=>\$bwa,
    "samtools=s"=>\$samtools,
    "annovar=s"=>\$annovar,
    "output=s"=>\$output,
    "help" => \$help,
);

die "$usage" if($help);
die "There is no RNA alignment file (rnabam) exists, please provide it!\n" if (!$rnabam);
die "There is no reference (hg19.fa) exists, please provide it!\n" if (!$reference);
die "There is no simpleRepeat file exists, please provide it!\n" if (!$simpleRepeat);
die "There is no alu file exists, please provide it!\n" if (!$alu);
die "There is no concat reference (hg19_mRNA.fa) exists, please provide it!\n" if (!$hg19_mRNA);
die "There is no bwa directory exists, please provide it!\n" if (!$bwa);
die "There is no samtools directory exists, please provide it!\n" if (!$samtools);
die "There is no annovar directory exists, please provide it!\n" if (!$annovar);
die "There is no output directory exists, please provide it!\n" if (!$output);

print "Begin time: ".`date`;

unless (-d $output) {`mkdir -p $output`;}

#Calculate mismatch ratio;
print "Calculating mismatch ratio...\n";
`$Bin/MismatchStat -i $rnabam -x 1000000 -u -o $output/rna_mismatch_stat.txt`;

#Calling varitation sites;
print "Variation sites detecting...\n";
`$Bin/MutDet -i $rnabam -r $reference  -q 50 -v $output/rna_mismatch_stat.txt -u -o $output/variation.sites.txt.gz`;

#Binomial,variation reads,editing frequency,strandbias,reads end,simple repeat and homopolymer filter;
print "Binomial,variation reads,editing frequency,strandbias,reads end,simple repeat and homopolymer filtering...\n";
`perl $Bin/binom.reads.fre.end.strand.mism.poly.rep.filter.pl --in $output/variation.sites.txt.gz --out $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.txt --simpleRepeat $simpleRepeat --reference $reference --endfre 0.9 --endp 0.05 --strandfre 0.9 --strandp 0.005 --homopolymer 5`;

#Merge the dbsnp files;
print "Reading and merging the snp database files...\n";
my @list;
open DBSNP, $snpdb or die $!;
while (<DBSNP>){
    chomp;
    push @list, $_;
}
close DBSNP;

open TEMP_SNP, ">$output/snp.temp.txt" or die $!;
foreach my $dbsnpfile (@list){
    if ($dbsnpfile =~ /\.gz$/){open SNP, "gzip -dc $dbsnpfile |" or die $!;}
    else {open SNP, $dbsnpfile or die $!;}
    while (<SNP>){
        chomp;
        if ($_ =~ /^#/) {next;}
        my @array = split /\s+/,$_;
        print TEMP_SNP $array[0]."\t".$array[1]."\n";
    }
    close SNP;
}

#dbSNP filter;
print "Removing the sites which in snp database...\n";
`perl $Bin/snp.filter.pl $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.txt $output/snp.temp.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt`;
`rm -f $output/snp.temp.txt`;

#Using the DNA bam file to filter SNP sites;
if (defined $dnabam){
    print "Using dna bam to filter SNP sites...\n";
    `$Bin/MismatchStat -i $dnabam -x 1000000 -o $output/dna_mismatch_stat.txt -u`;
    `perl $Bin/pileupSite.pl --bam $dnabam  --reg $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt --out $output/DNA.pileup --ref $reference --minMQ 20 --minBQ 20 --uniq`;
    `$Bin/pileupStat -i $output/DNA.pileup -q 20 -v $output/dna_mismatch_stat.txt -o $output/dna.mut.txt.gz`;
    `perl $Bin/filterDnaMut.pl $output/dna.mut.txt.gz > $output/dna.mut.filt.txt`;
    `perl $Bin/snp.filter.pl $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt $output/dna.mut.filt.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.dna.txt`;
}
if (defined $dnabam){`mv $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.dna.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt`;}
else {`mv $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt`;}

#Using bwa to remove false positive sites;
unless (-d "$output/bwa") {`mkdir -p $output/bwa`;}
print "Extract reads in which candidate RNA editing sites from the rna bam...\n";
`perl $Bin/extract_reads.pl $rnabam $output/bwa/mutation.read1.fq $output/bwa/mutation.read2.fq $output/bwa/mutation.read.fq $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt $output/readsvnum.txt`;
print "Realign the reads to the combined reference...\n";
`$bwa/bwa aln -o 1 -l 31 -k 2 -t 4 -L -i 15 $hg19_mRNA $output/bwa/mutation.read1.fq > $output/bwa/mutation.read1.sai`;
`$bwa/bwa aln -o 1 -l 31 -k 2 -t 4 -L -i 15 $hg19_mRNA $output/bwa/mutation.read2.fq > $output/bwa/mutation.read2.sai`;
`$bwa/bwa sampe -a 170*1.3 $hg19_mRNA $output/bwa/mutation.read1.sai $output/bwa/mutation.read2.sai $output/bwa/mutation.read1.fq $output/bwa/mutation.read2.fq |$samtools/samtools view -b -S - -t $hg19_mRNA.fai > $output/bwa/mutation.readP.bam`;
`$bwa/bwa aln -o 1 -e 50 -m 100000 -l 32 -k 2 -t 4 -L -i 15 $hg19_mRNA $output/bwa/mutation.read.fq > $output/bwa/mutation.read.sai`;
`$bwa/bwa samse $hg19_mRNA $output/bwa/mutation.read.sai $output/bwa/mutation.read.fq  | $samtools/samtools view -b -S -t $hg19_mRNA.fai - > $output/bwa/mutation.readS.bam`;
print "Using bwa mapping result to remove false positive sites...\n";
`perl $Bin/bwa.fre.filt.pl $output/bwa/mutation.readP.bam $output/bwa/mutation.readS.bam $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.txt $annovar/humandb_hg19/hg19_refGene.txt $output/readsvnum.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt`;

#RNA editing sites sorted by chr and position;
print "RNA editing sites sorting...\n";
`sort -k1.4n -k2n $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt > $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.sort.txt`;
`mv $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.sort.txt $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt`;

#Annovar annotation
print "Candidate RNA editing sites annotation...\n";
`perl $Bin/annovar.pl $output/binom.reads.fre.strand.end.simplerepeat.homopolymer.snp.bwa.txt $annovar $annovar/humandb_hg19/hg19_refGene.txt $alu $output/RNA_editing.sites.annotation.txt`;

#Splicing filter;
print "Splicing region filtering...\n";
`perl $Bin/alu.splicing.filter.pl $output/RNA_editing.sites.annotation.txt $output/RNA_editing.sites.annotation.splicing.filter.txt`;

#Basic information statistics;
print "Basic information counting...\n";
`perl $Bin/stat.pl $output/RNA_editing.sites.annotation.splicing.filter.txt $output/TypeDist.stat.txt $output/GenomeDist.stat.txt`;

#Delete temp files;
print "Deleting temp files...\n";
`rm -f $output/variation.sites.txt.gz`;
`rm -f $output/binom.reads.*`;
`rm -f $output/dna*`;
`rm -f $output/DNA*`;
`rm -rf $output/bwa`;
`rm -f $output/readsvnum*`;
`rm -f $output/rna_mismatch*`;
`mv $output/RNA_editing.sites.annotation.splicing.filter.txt $output/RNA_editing.sites.annotation.txt`;

print "End time: ".`date`;
