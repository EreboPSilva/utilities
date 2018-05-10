#!/usr/bin/perl

use strict;

my $infile = shift;
my $outfile = "fixed";
my $l = 60;

open (IN, $infile);
open (OUT, ">$outfile");
my $buffer = '';
my $line = '';
while (<IN>){
  unless (left($_, 1) eq '>'){
    chomp;
    $buffer .= $_;
    if (length($buffer) >= $l){
      $line = left($buffer, $l);
      $buffer = dleft($buffer, $l);
      print OUT "$line\n";
    }
  }
  else{
    print OUT join("\n", unpack( "(A$l)*", $buffer ))."\n";
    $buffer = '';
    print OUT $_;
    print $_;
  }
}
print "Done!\n";
close OUT;
close IN;


sub dright {
  # Deletes 'n' characters to the right of 'text'
  my $text   = shift;
  my $n      = shift;
  my $result = '';
  my $l      = length($text);
  if ( $n > 0 && $n < $l ) {
    $result = left( $text, $l - $n );
  }
  $result = $text if ( $n <= 0 );
  return $result;
}

sub dleft {
  # Deletes 'n' characters to the left of 'text'
  my $text   = shift;
  my $n      = shift;
  my $result = '';
  my $l      = length($text);
  if ( $n > 0 && $n < $l ) {
    $result = right( $text, $l - $n );
  }
  $result = $text if ( $n <= 0 );
  return $result;
}

sub left {
  # Returns 'n' characters to the left of 'string'
  my $string = shift;
  my $n      = shift;
  my $result = '';
  if ( $n > 0 ) {
    $result = substr( $string, 0, $n ) if ( $n < length($string) );
    $result = $string if ( $n >= length($string) );
  }
  return $result;
}

sub right {
  # Returns 'n' characters to the right of 'string'
  my $string = shift;
  my $n      = shift;
  my $result = '';
  if ( $n > 0 ) {
    $result = substr( $string, length($string) - $n, $n )
      if ( $n < length($string) );
    $result = $string if ( $n >= length($string) );
  }
  return $result;
}
