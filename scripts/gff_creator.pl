#!/usr/bin/perl -w
use strict;

#This script takes a list of .gff files, and prints STDOUT the new complete gff

my $fhlist = shift;
my $holder = slurp($fhlist);
my @allfiles = split /\n/, $holder;
my $IDnumber = 1;

for (my $p = 0; $p < scalar @allfiles; $p++){
my $gff;
my $fh = $allfiles[$p];

my $oldgff = slurp($fh);

$fh =~ s/\.gff//;
$fh =~ s/.+\///;

my @lines = split /\n/, $oldgff;
my @template = split /\t/, $lines[0];
my @last = split /\t/, $lines[-1];

my $c1 = $template[0];
my $c2 = $template[1];
my $c3;
my $from = $template[3];
my $to = $last[4];
my $c6 = $template[5];
my $c7 = $template[6];
my $c8 = '.';
my $gene = 'gene_' . sprintf("%04d",$IDnumber);
my $mRNA = 'mRNA_' . sprintf("%04d",$IDnumber);
my $exon = 'exon_' . sprintf("%04d",$IDnumber);
my $cds = 'cds_' . sprintf("%04d",$IDnumber);
my $name = $template[8];

$gff = join "\t", ($c1, $c2, "gene", $from, $to, $c6, $c7, $c8, "ID=$gene;Name=$name\n");
$gff .= join "\t", ($c1, $c2, "mRNA", $from, $to, $c6, $c7, $c8, "ID=$mRNA;Parent=$gene;Name=$name\n");

my $tmp;
for (my $j = 0; $j < scalar @lines; $j++){
    my @line = split /\t/, $lines[$j];
    $line[2] = 'CDS';
    $line[8] = "ID=$cds;Parent=$mRNA;Name=$name\n";
    $tmp .= join "\t", @line;
    $line[2] = 'exon';
    $line[7] = '.';
    $line[8] = "ID=$exon;Parent=$mRNA;Name=$name\n";
    $gff .= join "\t", @line;

}

$gff .= $tmp;

print $gff;

$IDnumber++;
}

print STDERR "Done!\n";

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
