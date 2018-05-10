#!/usr/bin/perl -w
use strict;

#IMPORTANT: Change the pattern for header-recognition in every case.
#	    Comment the rest, in order to preserve them.

my $fh = shift;
my $fasta = slurp($fh);
  chomp $fasta;
my $outf = $fh . ".fa";

$fasta =~ s/\n//g;
$fasta =~ s/\s/;/g;
$fasta =~ s/>/\n>/g;
$fasta =~ s/^\n//gs;

my $results;

my @allcontigs = split /\n/, $fasta;

for (my $i = 0; $i <= $#allcontigs; $i++){
  #my @contig = split /(>Contig\d+(?:;\[.+\])?)/, $allcontigs[$i];
  my @contig = split /(>..[^ATGCatgc]+)/, $allcontigs[$i];
  my $header = $contig[1];
  my $seq = $contig[2];
    $seq =~ s/(.{1,60})/$1\n/g;
  $results .= $header . "\n" . $seq . "\n";
}

open (OUT, '>', $outf) or die "Could not open $outf, $!";
print OUT $results;
close OUT;

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
