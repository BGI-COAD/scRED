#!usr/bin/perl -w
use strict;
die "perl $0 bamp bams tpv refgene bpv out" if (@ARGV != 6);
my ($in1,$in2,$in3,$in4,$in5,$out) = @ARGV;
open IN1, "samtools view $in1 |" or die $!;
open IN2, "samtools view $in2 |" or die $!;
open IN3, $in3 or die $!;
open IN4, $in4 or die $!;
open IN5, $in5 or die $!;
open OUT, ">$out" or die $!;
my @chr=();
for (my $i=1;$i<=22;$i++){push @chr,"chr$i";}
push @chr,"chrX";
push @chr,"chrY";
my %chr=();
foreach (@chr){$chr{$_}=1;}
my %hash;
my %hashall;
while (<IN4>){
    chomp;
    my @lines = split /\s+/,$_;
    if(exists $chr{$lines[2]}){
        my $tran = $lines[1];
        $hashall{$tran} = $_;
    }
}

sub map{
    my %hashrna;
    my ($lines,$pos) = @_;
    my @array = split /\s+/,$lines;
    my @start = split /,/,$array[9];
    my @end = split /,/,$array[10];
    my $len = @start;
    my $mrnapos = 1;
    for (my $i=1;$i<=$len;$i++){
        if ($array[3] eq "+"){
            for (my $j=$start[$i-1]+1;$j<=$end[$i-1];$j++){
                $hashrna{$mrnapos} = $j;
                $mrnapos += 1;
            }
        }
        if ($array[3] eq "-"){
            for (my $j=$end[-$i];$j>=$start[-$i]+1;$j--){
                $hashrna{$mrnapos} = $j;
                $mrnapos += 1;
            }
        }
    }
    return $hashrna{$pos};
}
while (<IN1>){
    chomp;
    if ($_ !~ /X0:i:1\t/){next;}
    my @line = split /\s+/,$_;
    my $name = $line[0];
    my $topchr;
    if ($name =~ /_(.+):/){$topchr = "chr".$1;}
    my @array = split /-/,$name;
    my $posv = $array[-4];
    my $readspos = $array[-3];
    my $pev = $array[-2];
    my $pos = $line[3];
    my $chr = $line[2];
    my $mapchr;
    my $rlen = length($line[9]);
    my $bin = unpack("B32",pack("N",$line[1]));
    my @bin = split //,$bin;
    @bin = reverse @bin;
    my $mappos;
    if (exists $hashall{$chr}){
        my @rnaline = split /\s+/,$hashall{$chr};
        $mapchr = $rnaline[2];
        if (exists $chr{$mapchr} and exists $chr{$topchr}){
            if ($bin[4] == 1){
                $readspos = $rlen - $readspos + 1;
                $mappos = $pos + $readspos -1;
                $pos = &map($hashall{$chr},$mappos);
            }
            else {
                $mappos = $pos + $readspos -1;
                $pos = &map($hashall{$chr},$mappos);
            }
        }
        else {next;}
    }
    else {
        $mapchr = $chr;
        if (exists $chr{$mapchr} and exists $chr{$topchr}){
            if ($bin[4] == 1){
                $readspos = $rlen - $readspos + 1;
                $pos = $pos + $readspos -1;
            }
            else {
                $pos = $pos + $readspos -1;
            }
        }
        else {next;}
    }
    my $pe;
    if ($bin[6]==1){$pe = 1;}
    if ($bin[7]==1){$pe = 2;}
    if(!$pos and $mapchr ne "*"){
        print STDERR "Warnning! Can not map to genome position: $name\t$topchr:$posv\t$mapchr:$mappos\.\n";
        next;
    }
    if ($mapchr eq "*"){
        print STDERR "Warnning! The reads can not map to reference!\n";
        next;
    }
    if ($pev == $pe and $posv ==$pos and $mapchr eq $topchr){
        if (exists $hash{$topchr."\t".$posv}){
            $hash{$topchr."\t".$posv} += 1;}
        else {
            $hash{$topchr."\t".$posv} = 1;
        }
    }
    else {
        if (exists $hash{$topchr."\t".$posv}){
            $hash{$topchr."\t".$posv} += 0;}
        else {
            $hash{$topchr."\t".$posv} = 0;
        }
    }
}
while (<IN2>){
      chomp;
      if ($_ !~ /X0:i:1\t/){next;}
      my @line = split /\s+/,$_;
      my $name = $line[0];
      my $topchr;
      if ($name =~ /_(.+):/){$topchr = "chr".$1;}
      my @array = split /-/,$name;
      my $posv = $array[-2];
      my $pos = $line[3];
      my $chr = $line[2];
      my $readspos = $array[-1];
      my $mapchr;
      my $rlen = length($line[9]);
      my $bin = unpack("B32",pack("N",$line[1]));
      my @bin = split //,$bin;
      @bin = reverse @bin;
      my $mappos;
      if (exists $hashall{$chr}){
          my @rnaline = split /\s+/,$hashall{$chr};
          $mapchr = $rnaline[2];
          if (exists $chr{$mapchr} and exists $chr{$topchr}){
              if ($bin[4] == 1){
                  $readspos = $rlen - $readspos + 1;
                  $mappos = $pos + $readspos -1;
                  $pos = &map($hashall{$chr},$mappos);
               }
               else {
                  $mappos = $pos + $readspos-1;
                  $pos = &map($hashall{$chr},$mappos);
               }
           }
           else {next;}
      }
      else {
          $mapchr = $chr;
          if (exists $chr{$mapchr} and exists $chr{$topchr}){
              if ($bin[4] == 1){
                  $readspos = $rlen - $readspos + 1;
                  $pos = $pos + $readspos -1;
              }
              else {
                  $pos = $pos + $readspos -1;
              }
          }
          else {next;}
      }
      if(!$pos and $mapchr ne "*"){
           print STDERR "Warnning! Can not map to genome position: $name\t$topchr:$posv\t$mapchr:$mappos\.\n";
           next;
      }
      if ($mapchr eq "*"){
           print STDERR "Warnning! The reads can not map to reference!\n";
            next;
      }
      if ($posv ==$pos and $mapchr eq $topchr){
           if (exists $hash{$topchr."\t".$posv}){
                $hash{$topchr."\t".$posv} += 1;}
           else {
                $hash{$topchr."\t".$posv} = 1;
           }
      }
      else {
           if (exists $hash{$topchr."\t".$posv}){
                $hash{$topchr."\t".$posv} += 0;}
           else {
                $hash{$topchr."\t".$posv} = 0;
           }
      }
}
my %hashtv;
while (<IN3>){
    chomp;
    my @array = split /\s+/,$_;
    my $chr = $array[0];
    my $pos = $array[1];
    my $vnum = $array[12];
    $hashtv{$chr."\t".$pos} = $_;
}
my %hashbv;
while (<IN5>){
    chomp;
    my @array = split /\t/,$_;
    my $chr = $array[0];
    my $pos = $array[1];
    my $vnum = $array[2];
    $hashbv{$chr."\t".$pos} = $vnum;
}
for my $key (keys %hash){
    if ($hashbv{$key} == 0) {next;}
    my $ratio = $hash{$key}/$hashbv{$key};
    if ($ratio >= 0.9){
        print OUT $hashtv{$key}."\t".$hashbv{$key}."\t".$hash{$key}."\t".$ratio."\n";
    }
}
close IN1;
close IN2;
close IN3;
close OUT;
