#! /usr/bin/perl -w 
use strict;
use Switch;

die "perl $0 <in.bam> <outpe1.fq> <outpe2.fq> <outse.fq> <var_region.txt> <var_support_stats>\nBAM files should have the flag NH:i:1\n" unless @ARGV==6;
my ($inbam, $outpe1, $outpe2, $outse, $invar, $outstats) = @ARGV;
my $insert ||= 500;
my %hash;
open OUTST, ">$outstats"or die $!;
open OUTPE1, ">$outpe1" or die $!;
open OUTPE2, ">$outpe2" or die $!;
open OUTSE, ">$outse" or die $!;

open VAR, $invar or die $!;
print "extract_reads is starting now!\n";
while (my $var = <VAR>)
{ #chr1	461390	8	A	4	4	0	4	0	1	3	C	4	2	2	3	1	0	4	2   0.428571	1	1
    chomp($var);
    my ($chr, $refPos, $refbase, $altbase) = (split (/\t/, $var))[0, 1, 3, 11];
    $refbase = substr($refbase,0,1);
    $altbase = substr($altbase,0,1);
    my $samPos = "$chr:$refPos-$refPos";
    my $st = $refPos - $insert;
    my $end = $refPos + $insert;
    my $regoin = "$chr:$st-$end";
    &reads ($refPos, $samPos, $refbase, $altbase, $regoin);
}

sub reads 
{
    my ($refPos, $samPos, $refBase, $altBase, $regoin) = @_;
    open IN, "samtools view  -h $inbam $samPos |" or die $!;
    my $total_reads = 0;
    my %temp;
    while (my $bamline = <IN>)
    {
        next if ($bamline =~ /^@/);
        chomp($bamline);
        my ($name, $flag, $chr, $alignPos, $MQ, $CIGAR, $peChr, $peAlignPos, $seq, $qual) = (split (/\t/,$bamline))[0, 1, 2, 3, 4, 5, 6, 7, 9, 10];
        if ($bamline =~ /\tNH:i:1/ and $MQ >= 50)
	{
	    my $readPos = &readloc_find( $refPos, $refBase, $altBase, $alignPos, $CIGAR, $seq, $qual);
#	    my ($readPos, $qbase, $basequa) = &readloc_find( $refPos, $refBase, $altBase, $alignPos, $CIGAR, $seq, $qual);   #$pos,$ref,$alt,$curpos,$cig,$seq
            if ( $readPos ne "down")
	    {
                my $qbase  = substr($seq, $readPos-1, 1);
                my $basequa = substr($qual, $readPos-1, 1);
                my $basequat = ord($basequa) - 33;
                if (($qbase ne $refBase) and ($qbase eq $altBase) )
		{
#			print "readPos:$readPos\tqbase:$qbase\trefBase:$refBase\taltBase:$altBase\n";
                    if(($flag & 16) == 16) #reads map to reverse strand
		    {
                        my @seq = split //,$seq;
                        my @qual = split //,$qual;
                        @seq = reverse @seq;
                        @qual = reverse @qual;
                        for (my $i=0; $i<@seq ; $i++) { $seq[$i] =~ tr/ACGTacgt/TGCAtgca/; }
                        $seq = join '',@seq;
                        $qual = join '',@qual;
                        $readPos = length($seq) + 1 - $readPos;
                     }
                     if (($flag & 8) == 8) #reads mate unmap
		     {
                         print OUTSE "\@".$name."_".$samPos."-".$readPos."\n".$seq."\n"."+"."\n".$qual."\n";
                         $total_reads ++;
#			print "$total_reads\tse\t$name\t$samPos $readPos $seq $qual\n";
                     }
	 	     if (($flag & 8) != 8) #reads mate mapped
#	 	     if ($peChr =~ /=/ or $peChr =~ /chr/)
		     {
                         my $pe = 0;
                         if (($flag & 64) == 64){$pe = 1;}  #first read in pe
                         if (($flag & 128) == 128){$pe = 2;} #second read in pe
                         my $pep = 0;
                         if (($flag & 64) == 64){$pep = 2;}
                         if (($flag & 128) == 128){$pep = 1;}
		 	 if ($pe != 0 && $pep != 0) {
			     $temp{$name}{$pep}{$samPos} = "\@".$name."_".$samPos."-\t$readPos\t-".$pe."\t".$seq."\n+\n".$qual."\n";
			     $total_reads ++;
#				print "$total_reads\tpe\t$name\t$samPos $readPos $seq \n";
			}
                    }
                }
            }
        }
    }
    close IN;
    my $outpos = (split /-/,$samPos)[1];
    my $outchr = (split /:/,$samPos)[0];
    print OUTST $outchr."\t".$outpos."\t".$total_reads."\n";

    open REG, "samtools view  $inbam $regoin|" or die $!;
    while(my $line = <REG>)
    {
    	chomp($line);
        if ($line !~ /^@/ && $line =~ /\tNH:i:1/)
        {
            my ($name, $flag, $readSeq, $readQua) = (split (/\t/,$line))[0, 1, 9, 10];
            my $pep=0;
            if (($flag & 64) == 64){$pep = 1;}
            if (($flag & 128) == 128){$pep = 2;}
            if (exists $temp{$name}{$pep})
            {
                if(($flag & 16) == 16)
                {
	            my @seq=split //,$readSeq;
                    my @qual=split //,$readQua;
                    @seq=reverse @seq;
                    @qual=reverse @qual;
                    for (my $i=0; $i<@seq ; $i++)
                    {
                    	$seq[$i]=~tr/ACGTacgt/TGCAtgca/;
                    }
                    $readSeq = join '',@seq;
                    $readQua = join '',@qual;
                }
                for my $key (keys %{$temp{$name}{$pep}})
                {
                    my @array = split /\t/,$temp{$name}{$pep}{$key};
                    my $pe = $array[2];
                    my $readPos = $array[1];
                    if (not exists  $array[0]) {next;}
                    if ($pep == 1)
		    {
                        print OUTPE1 "\@".$name."_".$key."-".$readPos.$pe."-1"."\n".$readSeq."\n"."+"."\n".$readQua."\n";
                        print OUTPE2 $array[0].$array[1].$array[2]."-2"."\n".$array[3];
                    }
                    if ($pep == 2)
                    {
                        print OUTPE2 "\@".$name."_".$key."-".$readPos.$pe."-2"."\n".$readSeq."\n"."+"."\n".$readQua."\n";
                        print OUTPE1 $array[0].$array[1].$array[2]."-1"."\n".$array[3];
                     }
                     $temp{$name}{$pep}{$key} =0;
                 }
             }
         }
    }
    close REG;

    for my $nam (keys %temp){
	for my $p (keys %{$temp{$nam}}){
	    for my $pos(keys %{$temp{$nam}{$p}}){
		$hash{$nam}{$p}{$pos} = $temp{$nam}{$p}{$pos} if ($temp{$nam}{$p}{$pos} ne "0");
	    }
	}
    }
    %temp = ();
}

open IN1, "samtools view  -h $inbam|" or die $!;
while (<IN1>){
    chomp;
    if ($_ !~ /^@/ && $_ =~ /\tNH:i:1/){
        my @lines = split /\t/,$_;
	my ($name, $flag) = ( @lines )[0, 1];
        my $pep=0;
        if (($flag & 64) == 64){$pep = 1;}
        if (($flag & 128) == 128){$pep = 2;}
        if (exists $hash{$name}{$pep}){
            if(($flag & 16) == 16){
                my @seq=split //,$lines[9];
                my @qual=split //,$lines[10];
                @seq=reverse @seq;
                @qual=reverse @qual;
                for (my $i=0; $i<@seq ; $i++){
                    $seq[$i]=~tr/ACGTacgt/TGCAtgca/;
                }
                $lines[9]=join '',@seq;
                $lines[10]=join '',@qual;
            }
            for my $key (keys %{$hash{$name}{$pep}}){
                my @array = split /\t/,$hash{$name}{$pep}{$key};
		#print "$name}{$pep}{$key\t$hash{$name}{$pep}{$key}\n";
                my $pe = $array[2];
                my $readPos = $array[1];
                if (not exists $array[0]) {next;}
                if ($pep == 1){
                    print OUTPE1 "\@".$name."_".$key."-".$readPos.$pe."-1"."\n".$lines[9]."\n"."+"."\n".$lines[10]."\n";
                    print OUTPE2 $array[0].$array[1].$array[2]."-2"."\n".$array[3];
                }
                if ($pep == 2){
                    print OUTPE2 "\@".$name."_".$key."-".$readPos.$pe."-2"."\n".$lines[9]."\n"."+"."\n".$lines[10]."\n";
                    print OUTPE1 $array[0].$array[1].$array[2]."-1"."\n".$array[3];
                }
            }
        }
    }
}
close OUTPE1;
close OUTPE2;
close OUTSE;
close IN1;

sub readloc_find 
{
    my ($pos, $ref, $alt, $curpos, $cig, $seq, $qual) = @_; #exonic SNP pos; reference base; alter base; reads mapping pos; CIGAR 12M1I23M; read sequence
#    print "var\t$pos, $ref, $alt, $curpos, $cig, $seq, $qual\n";
    my @arr1 = split(/\D/,$cig);
    my @arr2= split(/\d+/,$cig);
  
    my $readLen = length($seq);
    my $readloc = 1;
#    my $base="N";
#    my $Q="H";
    for( my $i=0; $i< @arr1; $i++)
    {
	switch ($arr2[$i+1])
    	{
	    case "M"
            {
                $curpos += $arr1[$i];
                $readloc += $arr1[$i];
                if ($curpos -1  >= $pos )
                {
                    $readloc -= ($curpos - $pos);
#                   $base = substr($seq, $readloc-1, 1) if ($readloc < $readLen);
                }
            }
            case "S" {$readloc += $arr1[$i];}
            case "D" {$curpos += $arr1[$i]; if ($curpos - 1 >= $pos ) {$readloc = "down"}}
            case "N" {$curpos += $arr1[$i]; if ($curpos - 1 >= $pos ) {$readloc = "down"}}
            case "I" {$readloc += $arr1[$i]; }
            case "H" {$readloc += 0;}
            else  { $readloc = "down"; print STDERR "Previous case not true"}
    	}
	last if (($curpos - 1 >= $pos) || ($readloc eq "down"));
    }
#    if ($readloc ne "down"){
#    	$base = substr($seq, $readloc-1, 1);
#	$Q =  substr($qual, $readloc-1, 1);
#    $readloc;
#    }
    return $readloc;
}
