#!/usr/bin/perl -w
use strict;

my $cases = shift;
my $lesserThan = shift;
  $lesserThan--;
my $mode = shift || 0;
my $result = [];

foreach my $j (1..$cases){
  mymain($result, $lesserThan);
  $lesserThan--;
}

if ($mode){
  foreach my $k (@{$result}){
    print $k . "\n";
  }
}
else{
  print join(', ', @{$result});
  print "\n";
}

sub mymain{
  my $result = shift;
  my $lesserThan = shift;
  my $n = myrnd($lesserThan);
    $n++;
  
  if ($result){
    foreach my $i (@{$result}){
      if ($i <= $n){
        $n++;
      }
    }
  }
  push @{$result}, $n;
  
  return $result;
}

sub myrnd{
  my $n = shift;
  return 0 unless ($n);
  my $done = 0;
  my $result = 0;
  while (!$done){
    $result = 0;
    my $t = $n;
    while ($t){
      $t = $t >> 1;
      my $r = rand();
      $result = $result << 1;
      $result++ if ($r >= 0.5);
    }
    $done = 1 if ($result <= $n);
  }
  return $result;
}