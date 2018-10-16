#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $fasta = slurp($fh);
  chomp $fasta;

my @contig = split />/, $fasta;
  shift @contig;

for (my $n = 0; $n < scalar @contig; $n++){
  my @chunck = split /\n/, $contig[$n];
  my $name = shift @chunck;
  #my $m = $n++;
  my $number = sprintf ("%06d", $n);
  my $outfile = 'LGEO' . $number . '.fsa';
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
