#!/usr/bin/perl -w
use strict;

my $file = shift;                             #name of the file with the sequence
my $position = shift;                         #position of interest
  $position--;                                  #adapts position to 'array-notation'
my $range = shift;                            #range of interest: absolute number or diferent to left-right
  my $min;
  my $max;
  if ($range =~ /-/){                         #adapts range to 'arry-notation'
    ($min, $max) = split /-/, $range;
    $min = $position - $min;
    $max += $position; 
  }
  else{
    $min = $position - $range;
    $max = $position + $range;
  }
  my @rang = ($min..$max);                    #range at last

my $myfile;
my $seq;
my $result;
my $warnmess;

open ($myfile, "<$file") or die ($!);

while (my $t = <$myfile>){                    #slurps the sequence without newlines or gene-names or stuff into a var
  chomp $t;
  if ($t =~ /^[A-Z]/){
    $seq .= $t;
  }
}

my @sepseq = split //, $seq;                  #turns the string with the sequence into an array with a residue per element

foreach my $j (@rang){                          #creates a string with all the elements acording to the previously created array
  if ($sepseq[$j] && $j >= '0'){                  #avoids 'circle text' and errors
    $result .= $sepseq[$j];
  }
  else{
    $warnmess = "Selected ranged outputs protein lenght.\n";    #informs of that
  }
}

if ($warnmess){                               #prints all the stuff needed
  print "$result\n\n$warnmess";
}
else{
  print $result;
}

#                                             #a cute whale
#         .-'
#   .'--./ /     _.---.
#    '-,  (__..-`       \
#       \          .     |
#        `,.__.   ,__.--/
#          '._/_.'___.-`
#
