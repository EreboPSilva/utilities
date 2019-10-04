#!/usr/bin/perl -w
use strict;

my $result;
my $mutg = 0;
my $mutt = 0;
my $muta = 0;
my $mutc = 0;
my $indet = 0;
my $correct;
my $review;
my $input = shift;
my $filename = shift;
#my $input = 'human_3_loc.out';
my $output1 = $filename . 'correct_reads.txt';
my $output2 = $filename . 'unclassified.txt';

open my $in, "<:encoding(utf8)", $input or die "$input: $!";
while (my $line = <$in>){
  chomp $line;

  my @elements = split /	/, $line;

  if ($elements[0] =~ m/CCCACGT/){
    $mutg = $mutg + $elements[1];
    $correct .= "$elements[0]\n";
  }
  elsif ($elements[0] =~ m/CCCACTT/){
    $mutt = $mutt + $elements[1];
  }
  elsif ($elements[0] =~ m/CCCACAT/){
    $muta = $muta + $elements[1];
  }
  elsif ($elements[0] =~ m/CCCACCT/){
    $mutc = $mutc + $elements[1];
  }
  else {
    $review .= "$elements[0]\t$elements[1]\n";
    $indet = $indet + $elements[1];
  }
}
close $in;

my $places = 4;
my $factor = 10**$places;
my $total = $muta + $mutc + $mutg + $mutt + $indet;
my $pcG = $mutg / $total * 100;
$pcG = int($pcG * $factor) / $factor;
my $pcT = $mutt / $total * 100;
$pcT = int($pcT * $factor) / $factor;
my $pcC = $mutc / $total * 100;
$pcC = int($pcC * $factor) / $factor;
my $pcA = $muta / $total * 100;
$pcA = int($pcA * $factor) / $factor;
my $pcI = $indet / $total * 100;
$pcI = int($pcI * $factor) / $factor;
#$result = $mutg / $total * 100;

open(my $out1, '>', $output1) or die "Could not open file '$output1' $!";
print $out1 $correct;
close $out1;

open(my $out2, '>', $output2) or die "Could not open file '$output2' $!";
print $out2 $review;
close $out2;

print "Lecturas mutadas a G: $mutg ($pcG%).\nLecturas mutadas a C: $mutc ($pcC%).\nLecturas mutadas a T: $mutt ($pcT%).\nLecturas mutadas a A: $muta ($pcA%).\nLecturas que no se adaptan al patrón: $indet ($pcI%).\n\n";
