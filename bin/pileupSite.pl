#!usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
#===============================================================================
#       AUTHOR:  Dongbing Liu, liudongbing@genomics.cn
#      CREATED:  2014/08/25
#===============================================================================
my $usage=<<USAGE;

Description
	Pileup the input bam file using the Bio::DB::Sam module
Parameter
	--bam	[STR] input bam file, required
	--reg	[str] position to pileup or a position list file, format as: chr pos, required
	--ref	[STR] input reference fasta file, required
	--out	[STR] the output file, defualt[STDOUT]
	--minMQ	[INT] minimal mapping quality, default[0]
	--maxMQ	[INT] maximal mapping quality, default[60]
	--minBQ	[INT] minimal base quality, default[0]
	--uniq	discard the non-unique mapping reads by NH(RNA mapped by tophat2) or H0(DNA mapped by bwa) tag, default[null].
	--help	give this information
Exmple 
	perl $0 --bam in.bam --reg site.txt --out out.pileup --ref hg19.fa
	perl $0 --bam in.bam --reg chr17:7571719 --out out.pileup --ref hg19.fa

USAGE

my ($bam,$reg,$out,$ref,$minMapQ,$minBaseQ,$maxMapQ,$uniq,$help);
GetOptions(
	"bam=s"=>\$bam,
	"reg=s"=>\$reg,
	"ref=s"=>\$ref,
	"out=s"=>\$out,
	"minMQ=i"=>\$minMapQ,
	"maxMQ=i"=>\$maxMapQ,
	"minBQ=i"=>\$minBaseQ,
	"uniq"=>\$uniq,
	"help" => \$help,
);
die "$usage" if(!$bam || !$reg || !$ref || $help);

if(!defined $minMapQ){$minMapQ=0;}
if(!defined $maxMapQ){$maxMapQ=60;}
if(!defined $minBaseQ){$minBaseQ=0;}

my $sam = Bio::DB::Sam->new(-bam  =>"$bam",
                            -fasta=>"$ref",
);

if(defined $out){open OUT,">$out" or die $!;}
if(-f $reg){
	my $process=0;
	open REG,"$reg" or die $!;
	while(<REG>){
		chomp;
		my @regLine=split /\s+/,$_;
		pileup($regLine[0],$regLine[1]);
		$process++;
		if($process%10000==0){print "Processing: $process\n";}
	}
	close REG;
}else{
	if($reg=~/(\w+):(\d+)/){
		pileup($1,$2);
	}
}
if(defined $out){close OUT;}

sub pileup{
    my ($chr,$pos) = @_;
	my @alignments = $sam->get_features_by_location(-seq_id => $chr,-start  => $pos, -end => $pos);
	my $effective_reads=0;
    my $bases="";
    my $quals="";
    my $mapQs="";
    my $qposs="";
    my $tposs="";
	my $misms="";
	my $qlens="";
	foreach my $b (@alignments) {
		my $qual=$b->qual;
		if($qual<$minMapQ || $qual>$maxMapQ){next;}
		if(defined $uniq){
			if($b->has_tag('NH')){
				my $nh_tag = $b->get_tag_values('NH');
				if($nh_tag=~/^\d+$/ && $nh_tag>1){next;} ## not unique mapping of tophat aligned dna reads
			}
			if($b->has_tag('H0')){
				my $h0_tag = $b->get_tag_values('H0');
				if($h0_tag=~/^\d+$/ && $h0_tag>1){next;} ## not unique mapping of bwa aligned dna reads
			}
            if($b->has_tag('IH')){ # added by Donby 2015/1/8
                my $nh_tag = $b->get_tag_values('IH');
                if($nh_tag=~/^\d+$/ && $nh_tag>1){next;} ## not unique mapping of tophat aligned dna reads
            }
		}

		my $cigar=$b->cigar_array;
		my $seq=$b->qseq;
		my @scores = $b->qscore;
		my $alignPos=$b->pos;
		my $qpos=&readPos($cigar,$alignPos,$pos);
		if($qpos==-1){next;} ## ref skip region
		$qpos--;
		my $qqual=$scores[$qpos];
		if($qqual<$minBaseQ){next;} ## base quality filter
		$effective_reads++;

		my $qbase  = substr($b->qseq,$qpos,1);
		$qbase=$b->reversed?lc($qbase):$qbase;
		$bases.=$qbase;
		$mapQs.=($effective_reads==1?$qual:",".$qual);
		$quals.=chr($qqual+33);
		my $tpos=$qpos+1;
        $tposs.=($effective_reads==1?$tpos:",".$tpos);
		my $qlen=$b->l_qseq;
		$qlens.=($effective_reads==1?$qlen:",".$qlen);

		my $nmism=0;
		if($b->has_tag('MD')){
			my $md_tag=$b->get_tag_values('MD');
			$md_tag=~s/\^[A-Z]+//g;
			$nmism= $md_tag=~tr/ACGT//;
		}
		$misms.=($effective_reads==1?$nmism:",".$nmism);

	}
	my $ref = $sam->segment($chr,$pos,$pos)->dna;
	$ref=uc($ref);
    if(defined $out){print OUT "$chr\t$pos\t$ref\t$effective_reads\t$bases\t$quals\t$mapQs\t$tposs\t$misms\t$qlens\n";}
	else{print "$chr\t$pos\t$ref\t$effective_reads\t$bases\t$quals\t$mapQs\t$tposs\t$misms\t$qlens\n";}
};

sub readPos{
	my ($cigar,$alignPos,$refPos)=@_;
	$alignPos++;
	my $readPos=0;
	my @cigar=@{$cigar};
	my $end;
	for (my $c=0; $c<@cigar; $c++) {
		my $char=$cigar[$c][0];
		my $cnum=$cigar[$c][1];
		if($char eq "S"){
			$readPos+=$cnum;
		}elsif($char eq "H"){
			next;
		}elsif($char eq "M"){
			$end=$alignPos+$cnum-1;
			if($end<$refPos){
				$readPos+=$cnum;
				$alignPos+=$cnum;
			}elsif($end>=$refPos){
				$readPos+=$refPos-$alignPos+1;
				last;
			}
		}elsif($char eq "I"){
			$readPos+=$cnum;
		}elsif($char eq "D"){
			$end=$alignPos+$cnum-1;
			if($end<$refPos){
				$alignPos+=$cnum;
			}else{
				$readPos=-1;
				last;
			}
		}elsif($char eq "N"){
			$end=$alignPos+$cnum-1;
			if($end<$refPos){
				$alignPos+=$cnum;
			}else{
				$readPos=-1;
				last;
			}
		}else{
			print STDERR "Cigar process error at $alignPos: $cigar\n";
		}
	}
	return ($readPos);
}
