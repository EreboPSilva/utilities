#!/usr/bin/perl -w
use strict;

#Este escript toma un genoma y escupe el mismo pero con el formato fasta bien hecho. 

my $fh = shift;
my $genome = slurp_esp($fh);
my $ofh = $fh . '_fixed';

open(my $OUT, '>', $ofh) or die "Could not open file '$ofh' $!";
print $OUT $genome;
close $OUT;

print "Done!\n";

sub slurp_esp{
  my $file = shift;
  my $tmp;
  my $result;

  open (my $fh, '<', $file) or die "cannot open file $file";

  while(my $line = <$fh>){
	if ($line =~ /^>Contig/){
		print "Processing: $line ...\n";
		$tmp =~ s/\n//g or print "Not possible to eliminate new lines in tmp, in $line";
		$tmp =~ s/(.{1,60})/$1\n/g or print "Not possible to edit lenght in tmp, in $line";
		$result .= $tmp or print "Not possible to include tmp in results, in $line";
		$tmp =~ s/.//g or print "Not possible to eliminate tmp, in $line";
		$result .= $line or print "Not possible to include contigs in results, in $line";
	}
	elsif($line =~ /^[ATGCN]/){
		$tmp .= $line; 
	}
	else{
		print "Weird line, check $line";
	}
  }
  
  close($fh);

  return $result;
}