#!/usr/bin/perl -w
use strict;

#Este escript toma Una lista de GOs separados por ';'.

my $inf = shift;
my $list = slurp($inf);
  chomp $list;
my %obj;
my $results;

my @tmp = split /;/, $list;

for my $i (0 .. $#tmp){
  $obj{$tmp[$i]}++;
}

my @ids = keys %obj;
for my $id (@ids){
  my $txt = $id . '; ' . $obj{$id} . "\n";
  $results .= $txt;
}

print $results;

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