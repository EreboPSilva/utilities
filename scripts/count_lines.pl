#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $li = slurp($fh);
  chomp $li;
my $lines;
my $result;

my @line = split /\n/, $li;
foreach my $m (@line){
    $lines->{$m}++;
}

my $n = 1;
foreach my $l (keys %{$lines}){
    $result .= $lines->{$l} . ' ' . $l . "\n";
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
