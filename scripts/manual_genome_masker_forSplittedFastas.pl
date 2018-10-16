#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $table = slurp($fh);
my @tab = split /\n/, $table;
my $coors = {};

for (my $i = 0; $i < scalar @tab; $i++){
    my @elements = split / {1,10}/, $tab[$i];
    $coors->{$elements[0]} = $elements[2];
}

foreach my $j (keys %{$coors}){
    my $fcontig = "fastas/$j.fsa.uncovered";
    my $contig = slurp($fcontig);
    my @lines = split /\n/, $contig; 
    my $result = shift @lines;
        $result .= "\n";
    my $joined = join '', @lines;
    my @letters = split //, $joined;

    my @initend = split /\.\./, $coors->{$j};
        $initend[0]++;
        $initend[1]++;
    my $diff = $initend[1] - $initend[0];
    my @xxx;
    for (my $k = 0; $k < $diff; $k++){
        push @xxx, 'X';
    }
    my @seq = splice @letters, $initend[0], $diff, @xxx;
    my $sequence = join '', @seq;
        $sequence =~ s/(.{1,60})/$1\n/g;
    $result .= $sequence;

    my $outf = 'fastas/' . $j . '.out';
    open (OUT, '>', $outf) or die "Couldn't open $outf";
    print OUT $result;
    close OUT;

    print "mv fastas/$j.fsa fastas/$j.fsa.uncovered\nmv fastas/$j.out fastas/$j.fsa\n";
}   

print "Done!\n";

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
