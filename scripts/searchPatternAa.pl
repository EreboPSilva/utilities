#!/usr/bin/perl -w
use strict;

my $pattern = shift;
my $fh = shift;
my $fasta = slurp($fh);
  chomp $fasta;

print $pattern;
print $fasta;

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
