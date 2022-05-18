#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use Getopt::Long;
my $usage=<<USAGE;

Description:
    This program is used to identify and filter RNA editing sites based on tophat alignment file(bam).  All of the programs were located in directory of bin when the software package was decompressed. The software package could only run on linux platform currently and the main program is RNA_editing_BGI.pl. 

Requirements:    
    Before running this program, you need to make sure that several pieces of software and database are installed on the system,
    BWA software, downloaded from:http://bio-bwa.sourceforge.net
    Samtools software, downloaded from:http://samtools.sourceforge.net
    Annovar, downloaded from:http://annovar.openbioinformatics.org
    Hg19 reference, downloaded from:http://genome.ucsc.edu
    Simple repeat region annotation file, downloaded from:http://genome.ucsc.edu
    Alu region annotation file,downloaded from:http://genome.ucsc.edu
    mRNA.fa, downloaded from:http://www.openbioinformatics.org
    dbSNP138, downloaded from:ftp://ftp.ncbi.nlm.nih.gov/snp.

Parameters:
    --rnabam       [STR] the RNA alignment file which is used to calling RNA editing sites.
    --reference    [STR] the genome fasta file(hg19.fa).
    --dbsnp        [STR] a file including SNP database files (including directory and file), eg ./temp/dbsnp138, each row indicate a SNP database.
    --simpleRepeat [STR] genome simple repeat region annotation file, should be bed format.
    --alu          [STR] genome alu region annotation file, should be bed format.
    --dnabam       [STR] the DNA aligment file, it could be used to remove SNP [optional].
    --hg19_mRNA    [STR] a combined reference (combination hg19.fa and mRNA.fa), and variation reads would be mapped to this reference using bwa.
    --bwa          [STR] the directory of bwa.
    --samtools     [STR] the directory of samtools.
    --annovar      [STR] the directory of annovar.
    --output       [STR] the directory of output.
    --help         [STR] show this information!

Example:
    Perl $0 --rnabam rnabam --reference hg19.fa --dbsnp dbsnp.txt --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --hg19_mRNA hg_mrna.fa --annovar annovar --output test

Notice:
    (1)The RNA and DNA alignment file (bam):
    The RNA bam file (tophat alignment) should be sorted, removed PCR duplicate, and we suggest you better to recalibrate base quality scores using gatk.
    The DNA bam file (BWA alignment) should also be sorted, removed PCR duplicate, and we suggest you better to recalibrate base quality scores and local indel realign using gatk.
    The rna and dna bam should be builded index, method: samtools index rna.bam.
    (2)Create the combined reference:
    Firstly, you should combine the hg19.fa and mRNA.fa; Secondly, build index for this modified reference using bwa, method: bwa index -a bwtsw hg19.fa.
    (3)The detail information of file format:
    The first two colums of SNP database file are chromsome and site, respectively.
    The format of Simple repeat region annotation file and Alu region annotation files is same, and the first three colums of this file are chromsome,start site and end site, respectively.
    (4)SNP filter:
    You must provide SNP database or dnabam to filter SNP at least, for example, you could only provide dbsnp138 if you do not have any other SNP database or dnabam. But we suggest that you should provide SNP datasets as many as possible in order to filter snp.

Outputs:
    When the program running completed, three files would be produced in the output directory.
    RNA_editing.sites.annotation.txt is RNA editing sites file, and the sites had been annotated by annovar software.
    TypeDist.stat.txt is the distribution of different kinds of RNA editing types, including content and number of each RNA editing type.
    GenomeDist.stat.txt is the distribution of RNA editing sites in genome, including content and number of RNA editing in each genomic region.

USAGE

my ($rnabam,$reference,$dbsnp,$simpleRepeat,$alu,$dnabam,$hg19_mRNA,$bwa,$samtools,$annovar,$output,$help);
GetOptions(
    "rnabam=s"=>\$rnabam,
    "reference=s"=>\$reference,
    "dbsnp=s"=>\$dbsnp,
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
die "$usage" if(!$rnabam || $help || !$output);

print "Begin time: ".`date`;

unless (-d $output) {`mkdir -p $output`;}
$reference ||= "/zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database/hg19.fa";
$dbsnp ||= "/zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database//dbsnp.list";
$simpleRepeat ||= "/zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database/hg19_simpleRepeat.reg.bed";
$alu ||= "/zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database/hg19.alu.bed";
$hg19_mRNA ||= "/zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database/hg_mrna/hg_mrna.fa";
$annovar ||= "/zfssz3/ST_CANCER/CGR/SHARE/RNA_editing_pipeline/hard_filter/database/annovar_lincRNA";
$bwa ||= "/zfssz3/ST_CANCER/CGR/USER/liudb2/bin/bwa-0.5.9";
$samtools ||= "/zfssz3/ST_CANCER/CGR/USER/liudb2/bin/samtools-0.1.18";

#Calculate mismatch ratio;
print "Calculating mismatch ratio...\n";
#`$Bin/MismatchStat -i $rnabam -x 1000000 -u -o $output/rna_mismatch_stat.txt`;
`$Bin/MismatchStat -i $rnabam  -u -o $output/rna_mismatch_stat.txt`;

#Calling varitation sites;
print "Variation sites detecting...\n";
`$Bin/MutDet -i $rnabam -r $reference  -q 50 -v $output/rna_mismatch_stat.txt -u -o $output/var.sites.txt.gz`;

#Binomial,variation reads,editing frequency,strandbias,reads end,simple repeat and homopolymer filter;
print "Binomial,variation reads,editing frequency,strandbias,reads end,simple repeat and homopolymer filtering...\n";
`perl $Bin/binom.reads.fre.end.strand.mism.poly.rep.filter.pl --in $output/var.sites.txt.gz --out $output/var.basic.txt --simpleRepeat $simpleRepeat --reference $reference --endfre 0.9 --endp 0.05 --strandfre 0.9 --strandp 0.005 --homopolymer 5`;

#Merge the dbsnp files;
print "Reading and merging the snp database files...\n";
my @list;
open DBSNP, $dbsnp or die $!;
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
`perl $Bin/snp.filter.pl $output/var.basic.txt $output/snp.temp.txt $output/var.basic.dbsnp.txt`;
`rm -f $output/snp.temp.txt`;

#Using the DNA bam file to filter SNP sites;
if (defined $dnabam){
    print "Using dna bam to filter SNP sites...\n";
    `$Bin/MismatchStat -i $dnabam -x 1000000 -o $output/dna_mismatch_stat.txt -u`;
    `perl $Bin/pileupSite.pl --bam $dnabam  --reg $output/var.basic.dbsnp.txt --out $output/DNA.pileup --ref $reference --minMQ 20 --minBQ 20 --uniq`;
    `$Bin/pileupStat -i $output/DNA.pileup -q 20 -v $output/dna_mismatch_stat.txt -o $output/dna.mut.txt.gz`;
    `perl $Bin/filterDnaMut.pl $output/dna.mut.txt.gz > $output/dna.mut.filt.txt`;
    `perl $Bin/snp.filter.pl $output/var.basic.dbsnp.txt $output/dna.mut.filt.txt $output/var.basic.dbsnp.dna.txt`;
}
if (defined $dnabam){`mv $output/var.basic.dbsnp.dna.txt $output/var.basic.snp.txt`;}
else {`mv $output/var.basic.dbsnp.txt $output/var.basic.snp.txt`;}

#Using bwa to remove false positive sites;
unless (-d "$output/bwa") {`mkdir -p $output/bwa`;}
print "Extract reads in which candidate RNA editing sites from the rna bam...\n";
`perl $Bin/extract_reads.pl $rnabam $output/bwa/mutation.read1.fq $output/bwa/mutation.read2.fq $output/bwa/mutation.read.fq $output/var.basic.snp.txt $output/readsvnum.txt`;
print "Realign the reads to the combined reference...\n";
`$bwa/bwa aln -o 1 -l 31 -k 2 -t 4 -L -i 15 $hg19_mRNA $output/bwa/mutation.read1.fq > $output/bwa/mutation.read1.sai`;
`$bwa/bwa aln -o 1 -l 31 -k 2 -t 4 -L -i 15 $hg19_mRNA $output/bwa/mutation.read2.fq > $output/bwa/mutation.read2.sai`;
`$bwa/bwa sampe -a 170*1.3 $hg19_mRNA $output/bwa/mutation.read1.sai $output/bwa/mutation.read2.sai $output/bwa/mutation.read1.fq $output/bwa/mutation.read2.fq |$samtools/samtools view -b -S - -t $hg19_mRNA.fai > $output/bwa/mutation.readP.bam`;
`$bwa/bwa aln -o 1 -e 50 -m 100000 -l 32 -k 2 -t 4 -L -i 15 $hg19_mRNA $output/bwa/mutation.read.fq > $output/bwa/mutation.read.sai`;
`$bwa/bwa samse $hg19_mRNA $output/bwa/mutation.read.sai $output/bwa/mutation.read.fq  | $samtools/samtools view -b -S -t $hg19_mRNA.fai - > $output/bwa/mutation.readS.bam`;
print "Using bwa mapping result to remove false positive sites...\n";
`perl $Bin/bwa.fre.filt.pl $output/bwa/mutation.readP.bam $output/bwa/mutation.readS.bam $output/var.basic.snp.txt $annovar/humandb_hg19/hg19_refGene.txt $output/readsvnum.txt $output/var.basic.snp.bwa.txt`;
`msort -k m1[4-5] -k n2  $output/var.basic.snp.bwa.txt > $output/var.basic.snp.bwa.sort.txt`;
`mv $output/var.basic.snp.bwa.sort.txt $output/var.basic.snp.bwa.txt`;

#Annovar annotation
print "Candidate RNA editing sites annotation...\n";
`perl $Bin/annovar.pl $output/var.basic.snp.bwa.txt $annovar $annovar/humandb_hg19/hg19_refGene.txt $alu $output/RNA_editing.sites.annotation.txt`;

#Splicing filter;
print "Splicing region filtering...\n";
`perl $Bin/alu.splicing.filter.pl $output/RNA_editing.sites.annotation.txt $output/RNA_editing.sites.annotation.splicing.filter.txt`;

#Basic information statistics;
print "Basic information counting...\n";
`perl $Bin/stat.pl $output/RNA_editing.sites.annotation.splicing.filter.txt $output/TypeDist.stat.txt $output/GenomeDist.stat.txt`;

#Delete temp files;
#print "Deleting temp files...\n";
#`rm -f $output/variation.sites.txt.gz`;
#`rm -f $output/binom.reads.*`;
#`rm -f $output/dna*`;
#`rm -f $output/DNA*`;
#`rm -rf $output/bwa`;
#`rm -f $output/readsvnum*`;
#`rm -f $output/rna_mismatch*`;
`mv $output/RNA_editing.sites.annotation.splicing.filter.txt $output/RNA_editing.sites.annotation.txt`;

print "End time: ".`date`;
