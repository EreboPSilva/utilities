#!/usr/bin/perl -w
use strict;

my $fh = shift;
my @fasta = slurp_array($fh);
	chomp @fasta;

my $header = shift @fasta;
	$header =~ s/>/>Geo_/;
my $seq = join '', @fasta;
	$seq =~ s/(.{1,60})/$1\n/g;
my $result = $header . "\n" . $seq;

print $result;

sub slurp_array{
  my $file = shift;
	my @result;
	
  open (IN, $file) or warn $file;
  foreach my $l (<IN>){
	  push @result, $l;
	}
  close IN;
  
  return @result;
}

sub slurp{
  my $file = shift;

  open (IN, $file) or warn $file;
  my $result = <IN>;
  close IN;
  
  return $result;
}