#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $seqs = slurp($fh);
  chomp $seqs;
my $sequences;
my $result;

my @seqline = split /\n/, $seqs;
foreach my $m (@seqline){
  $sequences->{$m}++;
}

my $n = 1;
foreach my $l (keys %{$sequences}){
  my $header = '>' . $sequences->{$l} . '.' . $n;
  $result .= $header . "\n" . $l . "\n";
  $n++;
}

print $result;

sub slurp{
  my $file = shift;
  my $result;

  open (my $sfh, '<', $file) or die "cannot open file $file";
  {
    local $/;
    $result = <$sfh>;
  }
  close($sfh);

  return $result;
}
