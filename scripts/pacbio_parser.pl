#!/usr/bin/env perl -w
use strict;
use Carp;

my $fhgfs = shift;
my $gapfillstatus = slurp($fhgfs);
my @liststatus = split /\n/, $gapfillstatus;
my @seqstoblast;
my $result;
my $log;
my $count = 0; 

print STDERR "Parsing " . scalar @liststatus . " entries and retrieving sequences...\n";
for (my $j = 0; $j < scalar @liststatus; $j++){
    my @element = split /\t/, $liststatus[$j];
    
    if (    $element[1] eq 'minreadfail' || 
            $element[1] eq 'filled' || 
            $element[1] eq 'nofillmetrics'){
        $count++;
    }
    elsif ( $element[1] eq 'singleextend' || 
            $element[1] eq 'doubleextend' || 
            $element[1] eq 'overfilled' ){
       my $seq = seqfinder($element[0]);
       push @seqstoblast, $seq;
    }
    else{
        print STDERR "Unknown status: $element[1].\n";
    }
}

print STDERR "Performing " . scalar @seqstoblast . " BLASTNs and extracting coordinates...\n\n";
for (my $k = 0; $k < scalar @seqstoblast; $k++){
    my $blast = `bash blaster.sh $seqstoblast[$k]`;
    if ($blast){
        my @q = split /\n/, $blast;
        my @w = split /\t/, $q[0];
        if ($q[1]){
            my @e = split /\t/, $q[1];
            my @coors = ($w[5], $w[6], $e[5], $e[6]);
                sort @coors;
        
            if ($w[3] == $w[10] && $w[4] eq $e[4] && $e[3] == $e[10]){
                $result .= $w[4] . "\t" . $coors[0] . "\t" . $coors[3] . "\n";
            }
            else{
                if ($w[1] == 1 && $w[3] == $w[10]){
                    $result .= $w[4] . "\t" . $w[5] . "\t" . $w[6] . "\n";
                }
                else{
                    print STDERR "There is not completness in first. Check: " . "\n" . $blast . "\n";
                }
            }
        }
        else{
            if ($w[3] == $w[10]){
                $result .= $w[4] . "\t" . $w[5] . "\t" . $w[6] . "\n";
            }
            else{
                print STDERR "There is not completness. Check: " . "\n" . $blast . "\n";
            }
        }
        $log .= $k . "\t" . $blast;
    }
    else{
        print STDERR "No blast result for $seqstoblast[$k]\n";
    }
}

print $result;

print STDERR "\n$count elements ignored due lack of assembly evidence.\n";

open(OUT, '>', 'pacbio_parser.log') or die "cannot open file OUT";
print OUT $log;
close OUT;

print STDERR "Done!\n";

sub slurp{
    my $file = shift;
    my $result;
    open (my $sfh, '<', $file) or die "cannot open file $file";
    {
        local $/;
        $result = <$sfh>;
    }
    close($sfh);
    return $result;
}

sub seqfinder{
    my $q = shift;
    my $fh = "assembly/$q/fillingMetrics.json";
    my $json = slurp($fh);

    my $seq = $1 if $json =~ /"extendSeq.": "([ATGC]+)",/;
        chomp $seq;

    return $seq;
}
