#!/usr/bin/perl -w
use strict;

#Este escript toma una serie de genomas, indicados mediante sus rutas absolutas separadas por comas, y debuelve una tabla indicando el contig (detonado solo con un numero) y la longitud, de forma que se pueda hacer un grafico de lineas con el.

my $fh = shift;
my @paths = split /,/, $fh;
my $n_genomes = scalar @paths;

my @genomes;
my @results;

for (my $i = 0; $i < $n_genomes; $i++){
  my $fasta = slurp($paths[$i]);
  push @genomes, $fasta;
}

for (my $j = 0; $j < $n_genomes; $j++){
  my @lengths = stracter($genomes[$j]);
  splice @results, 1, 0, \@lengths;
}

my @maxish;
for (my $d = 0; $d < $n_genomes; $d++){
  my $temp = scalar @{$results[$d]};
  push @maxish, $temp;
}
my @max = sort { $a <=> $b } @maxish;
my $last = pop @max;
$last++;

warn join ', ', @maxish;
warn join ', ', @max;

my $result = '# ' . $fh . "\n";

for (my $s = 1; $s <= $last; $s++){
  $result .= $s;
  for (my $q = 0; $q < $n_genomes; $q++){
    $result .= '; ' . $results[$q][$s];
  }
  $result .= "\n";
}

print $result;

sub stracter{
  my @result;
  my $genome = shift;
    chomp $genome;
  my @contigs = split />/, $genome;
    shift @contigs;
  my @joined_contigs;
  my @unsorted_result;
  
  foreach my $m (@contigs){
    my @seqs = split /\n/, $m;
    shift @seqs;
    my $tmp = join '', @seqs;
    push @joined_contigs, $tmp;
  }
  
  foreach my $g (@joined_contigs){
    my $length = length $g;
    push @unsorted_result, $length;
  }
  
  @result = sort { $b <=> $a } @unsorted_result;
  return @result;
}

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
}0
