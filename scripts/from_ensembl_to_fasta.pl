
#!/usr/bin/perl -w
use strict;

my $fh = shift;
my $header = $fh;
	$header =~ s/.fa//g;
	$header =~ s/^/>/;
my $fasta = slurp($fh);
	$fasta =~ s/^>.+\n//mg;
	$fasta =~ s/\n//g;
	$fasta =~ s/(.{1,60})/$1\n/g;

print $header;
print "\n";
print $fasta;

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
  my $result;

  open (IN, $file) or warn "$! [$file]";
  while (my $line = <IN>){
	chomp $line;
	$result .= $line . "\n";
  }
  close IN;
  
  return $result;
}