#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $fasta = slurp($fh);
    chomp $fasta;

my @contig = split />/, $fasta;
    shift @contig;

for (my $n = 0; $n < scalar @contig; $n++){
    my @chunck = split /\n/, $contig[$n];
    my $header = shift @chunck;
    my @words = split /\|/, $header;
    my $seq = join '', @chunck;
    my @sequences = split //, $seq;
    my $length = scalar @sequences;
    my $number = $n++;

    my $entrynumber = "CABI" . $number;
    my $organism = "Chelonoidis abingdonii";
    my $envsam = 'no';
    my $gene;
    my $product;
    my $translatable = ;
    my $cds5 = '1';
    my $cds3 = $length++;
    my $partial5  = 'yes';
    my $partial3 = 'yes'; ###
    my $transltable = '1';
    my $finalsequence = join '', @sequences;

    if ($words[4] eq 'unknown'){
        $gene = "Unknown";
        $product = "Protein of unknown function";
    }
    else {
        $gene = $words[4];
        for (my $m = 5; $m < scalar @words; $m++){
            $product .= $words[$m] . ' ' unless $words[$m] =~ m/[\(\)]/;
        }
        if (!$product){
            $product = $gene . " ";
        }
        $product =~ s/ $//;
    }

    if ($sequences[0] eq 'M'){
        $partial5 = 'no';
    }
    print $entrynumber . "\t" . $organism . "\t" . $gene . "\t" . $product . "\t" . $cds5 . "\t" . $cds3 . "\t" . $partial5 . "\t" . $partial3 . "\t" . $transltable . "\t" . $finalsequence . "\n";
}

print STDERR "Done!";

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
