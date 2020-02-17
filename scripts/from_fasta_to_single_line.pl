#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $fasta = slurp($fh);
  chomp $fasta;
my $result;

my @contig = split />/, $fasta;
    shift @contig;

for (my $n = 0; $n < scalar @contig; $n++){
    my @chunck = split /\n/, $contig[$n];
    my $name = shift @chunck;
    my $seq = join '', @chunck;
    $result .= '>' . $name . "\n" . $seq . "\n";
}

print $result;
warn "Done!";

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
