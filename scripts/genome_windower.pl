#!/usr/bin/perl -w
use strict;

#Este escript toma un genoma y escupe el mismo pero con "contigs" de Xkb de largo. 

my $fh = shift;
my $genome = slurp($fh);
  chomp $genome;
my $l = shift;
  my $length = $l * 1000;
my $result;

(my $bulk = $genome) =~ s/>Contig//g;
  $bulk =~ s/[0-9]//g;
  $bulk =~ s/\n//g;
  
(my $separated = $bulk) =~ s/(.{1,$length})/$1\n/g;
  my @windowed = split /\n/, $separated;
  shift @windowed;
my $n = scalar @windowed;

for (my $k = 0; $k < $n; $k++){
  $result .= '>' . $k . "\n";
  (my $tmp = $windowed[$k]) =~ s/(.{1,60})/$1\n/g;
  $result .= $tmp . "\n";
}

my @filename = split /\./, $fh;
my $outf = $filename[0] . '_w' . $l . 'kb.' . $filename[1];
open(my $OUT, '>', $outf) or die "Could not open file '$outf' $!";
print $OUT $result;
close $OUT;

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