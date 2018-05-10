#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $fasta = slurp($fh);
  chomp $fasta;

my @contig = split />/, $fasta;
  shift @contig;

foreach my $n (@contig){
  my @chunck = split /\n/, $n;
  my $name = shift @chunck;
  my $outfile = $name . '.fa';
  my $seq = join "\n", @chunck;
  my $text = '>' . $name . "\n" . $seq . "\n";
  
  open(my $foh, '>', $outfile) or die "Could not open file '$outfile' $!";
    print $foh $text;
  close $foh;
}

print "Done!";

sub slurp{
  my $file = shift;
  my $result;

  open (my $fh, '<', $file) or die "cannot open file $file";
  {
    local $/;
    $result = <$fh>;
  }
  close($fh);

  return $result;
}