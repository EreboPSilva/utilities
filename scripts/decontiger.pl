#!/usr/bin/perl -w
use strict;
use bytes;
my $N = shift || 20;
my $genome = shift || '/home/jmgps/testing/split_george_genome/cnig_sspace_lgeorge.fa';
#my $text = getfile($genome);
#my @contigs = split />/, $text;

#foreach my $k (@contigs){
warn("Checking Size...\n");
my $size = -s $genome;
my $max_size = int($size / $N);
my $fsize = 0;
my $fcounter = 0;
open (IN, $genome);
my $current_file = fname($fcounter);
open (OUT, ">$current_file");
while (<IN>){
  #warn "In the loop...";
  if (/^>(.+)/ && $fsize > $max_size){
    $fsize = 0;
    close OUT;
    $fcounter++;
    $current_file = fname($fcounter);
    open(OUT, ">$current_file");
    warn("Starting file $fcounter with contig $1");
  }
  print OUT $_;
  $fsize += length($_);
}
close OUT;
close IN;

sub fname{
  my $n = shift;
  my $file = "cnig_split_" . "$n" . ".fa";
  return $file;
}
sub getfile{
  warn "I'm here...";
  my $filename = shift;
  my $content;
  
  open IN, "< $filename" or die "Can't Open $filename: $!";
  while (my $j = <IN>){
    $content .= $j;
  }
  #{ local $/ = undef;
  #  $content = <IN>;
  #}
  close IN;
  
  return $content;  
}
