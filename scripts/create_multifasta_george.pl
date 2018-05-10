#!/usr/bin/perl -w
use strict;
use LWP;

my $ua = LWP::UserAgent->new;
$ua->agent('Mozilla/5.0 (Windows NT 10.0; WOW64; rv:50.0) Gecko/20100101 Firefox/50.0');
my $url = 'https://www.ncbi.nlm.nih.gov/protein/XP_006124860.1?report=fasta&log$=seqview&format=text';
my $response = $ua->get($url);
print($response->{'_content'});
#print join("\n", keys %$prev_links);
print "\nDone!\n";

