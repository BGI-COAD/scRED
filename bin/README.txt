Description:
    This software package is used to identify and filter RNA editing sites based on tophat alignment file(bam). All of the programs were located in directory of bin when the software package was decompressed. The software package could only run on linux platform currently and the main program is RNA_editing.pl.

Requirements:    
    Before running this program, you need to make sure that several pieces of software and database are installed on the system,
    BWA software, downloaded from:http://bio-bwa.sourceforge.net
    Samtools software, downloaded from:http://samtools.sourceforge.net
    Annovar, downloaded from:http://annovar.openbioinformatics.org
    Hg19 reference, downloaded from:http://genome.ucsc.edu
    Simple repeat region annotation file, downloaded from:http://genome.ucsc.edu
    Alu region annotation file,downloaded from:https://github.com
    mRNA.fa, downloaded from:http://www.openbioinformatics.org
    dbSNP138, downloaded from:ftp://ftp.ncbi.nlm.nih.gov/snp.

Parameters:
    --rnabam       [STR] the RNA alignment file which is used to calling RNA editing sites (required).
    --reference    [STR] the genome fasta file(hg19.fa) (required).
    --snpdb        [STR] a file including SNP database files (including directory and file), eg ./temp/dbsnp138, each row indicate a SNP database (required).
    --simpleRepeat [STR] genome simple repeat region annotation file, should be bed format (required).
    --alu          [STR] genome alu region annotation file, should be bed format (required).
    --dnabam       [STR] the DNA aligment file, it could be used to remove SNP [optional].
    --hg19_mRNA    [STR] a combined reference (combination hg19.fa and mRNA.fa), and variation reads would be mapped to this reference using bwa (required).
    --bwa          [STR] the directory of bwa (required).
    --samtools     [STR] the directory of samtools (required).
    --annovar      [STR] the directory of annovar (required).
    --output       [STR] the output directory (required).
    --help         [STR] show this information!

Example:
    Perl RNA_editing.pl --rnabam rnabam --reference hg19.fa --dbsnp dbsnp.txt --simpleRepeat hg19_simpleRepeat.reg.bed --alu hg19.alu.bed --hg19_mRNA hg_mrna.fa --annovar annovar --output test
    Perl RNA_editing.pl --help 

Notice:
    (1)The RNA and DNA alignment file (bam):
    The RNA bam file (tophat alignment) should be sorted, removed PCR duplicate, and we suggest you had better to perform recalibrate base quality scores using gatk.
    The DNA bam file (BWA alignment) should also be sorted, removed PCR duplicate, and we suggest you had better to perform recalibrate base quality scores and local indel realign using gatk.
    The rna and dna bam should be builded index, method: samtools index rna.bam.
    (2)Create the concat reference:
    Firstly, you should merge the hg19.fa and mRNA.fa; Secondly, build index for this concat reference using bwa, method: bwa index -a bwtsw hg19_mRNA.fa.
    (3)The detail information of file format:
    The first two colums of SNP database file are chromsome and position, respectively.
    The format of Simple repeat region annotation file and Alu region annotation files is bed format, and the first three colums of this file are chromsome,start position and end position, respectively.
    (4)SNP filter:
    You must provide SNP database or dnabam to filter SNP at least, for example, you could only provide dbsnp138 if you do not have any other SNP database or dnabam.

Outputs:
    When the program running completed, three files would be produced in the output directory.
    RNA_editing.sites.annotation.txt is RNA editing sites file, and the sites had been annotated by annovar software and alu region annotation files.
    TypeDist.stat.txt is the distribution of different kinds of RNA editing types, including content and number of each RNA editing type.
    GenomeDist.stat.txt is the distribution of RNA editing sites in the genome, including content and number of RNA editing in each genomic region.
