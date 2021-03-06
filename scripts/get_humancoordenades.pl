#!/usr/bin/perl -w
use strict;

my $name = shift;
my $fh = '/mnt/nas1/vqf/lgeorge/PAML/protein_databases/some_gaps/' . $name . '.cdna.fa.paml.best.phy.nsites.out';
my $pamloutput = slurp($fh);
my ($seqsh) = $pamloutput =~ m/(Hum_$name +.+)\n/;
  my @codonsh = split /\s+/, $seqsh;
  my @aminoacidsh = translateArray(@codonsh);
    unshift @aminoacidsh, $codonsh[0];
my ($seqsg) = $pamloutput =~ m/(Hum_$name +.+)\n/;
  my @codonsg = split /\s+/, $seqsg;
  my @aminoacidsg = translateArray(@codonsg);
#    unshift @aminoacidsg, $codonsg[0];
my ($residues) = $pamloutput =~ m/Bayes Empirical Bayes \(BEB\) analysis \(Yang, Wong & Nielsen 2005\. Mol\. Biol\. Evol\. 22:1107-1118\)\nPositive sites for foreground lineages Prob\(w>1\):\s+(.+)\n\n\nThe grid/s;
  $residues =~ s/\n +/\n/g;

my @variants = split /\n/, $residues;
foreach my $i (@variants){
  my @elements = split / /, $i;
  my $pos = $elements[0];
  my $res = $elements[1];
  
  shift @codonsg;
  unshift @codonsg, $pos;
  my $coor = getCoor(@codonsg);

  print $coor . "\t" . $name . '(' . $res . $pos . ')' . "\n";
#  print $coor . "\n"
}

sub getCoor{
  my @seq = @_;
  my $position = $seq[0];
  my $diff = 10;
    my $lastpos = $position + $diff;
	my $firstpos = $position - $diff;
  my @rang = ($firstpos..$lastpos);
  my $naseq;
  my $results;
  
  foreach my $k (@rang){
    if ($seq[$k]){
      $naseq .= $seq[$k];
    }
  }
  
  $naseq =~ s/---/nnn/g;
  
  system "/home/vqf/software/blatMeHuman.sh $naseq > /mnt/nas1/jmgps/tmp/get_coor.tmp";
  
  my $blatfile = '/mnt/nas1/jmgps/tmp/get_coor.tmp';
  my $blatresult = slurp($blatfile);
  
  my ($contig) = $blatresult =~ m/>(chr\d+)/;
  my ($coord) = $blatresult =~ m/Sbjct: (\d+)/;
  my ($direction) = $blatresult =~ m/Strand = Plus \/ (.+)\n/;
  my ($start) = $blatresult =~ m/Query: (\d*)/;
  my $number_to_operate = 30;
    unless ($start == 1){
	  $number_to_operate -= $start;
	}
  my $symbol;
  
  if ($direction eq 'Minus'){
    $coord -= $number_to_operate;
	$symbol = '-';
  }
  elsif ($direction eq 'Plus'){
    $coord += $number_to_operate;
	$symbol = '+';
  }
  
  my $base_change;
  my ($blat_bases_string) = $blatresult =~ m/  (.[agct-]+) /;
  my @blat_bases = split //, $blat_bases_string;
  my $index_blat = 30;
    $index_blat -= $start;
  my $original_base = $blat_bases[$index_blat];
  my $new_base;
  
  if ($original_base eq 't'){
    $new_base = 'c';
  }
  else{
    $new_base = 't';
  }
  
  $base_change = $original_base . "/" . $new_base;

#  $results = $contig . "\t" . $coord . "\t" . $coord . "\t" . $base_change . "\t" . $symbol;
  $results = $contig . $coord;

#  return $blatresult;
#  return $direction;
#  return $naseq;
  return $results;
}

sub translateArray{
  my @codons = @_;
    shift @codons;
  my(%g)=(
    'TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S',
    'TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L',
    'TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_',
    'TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W',
    'CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L',
    'CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P',
    'CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q',
    'CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R',
    'ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M',
    'ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T',
    'AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K',
    'AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R',
    'GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V',
    'GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A',
    'GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E',
    'GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G',
	'---'=>'X'
  );
  my @results;
  
  foreach my $i (@codons){
    if (exists $g{$i}){
	  push @results, $g{$i};
	}
	else{
	  print STDERR "Bad codon \"$i\"!!\n";
	  exit;
	}
  }
  return @results;
}

sub slurp{
  my $file = shift;
  my $result;

  open (IN, $file) or warn "$! [$file]";
  while (my $line = <IN>){
       chomp $line;
       $result .= $line . "\n";
  }
  close IN;

  return $result;
}