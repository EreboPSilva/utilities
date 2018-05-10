#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

		#tie %sequences,'Bio::DB::Fasta','/path/to/fasta/files';
		#print $sequences{'CHROMOSOME_I:1,20000'};

my $gff_file = shift || 'geosplitted.gff';
my $genome = shift || 'geosplitted.fa';

print "Tying genome...\n";
my %db = ();
my $tobj = tie %db, 'Bio::DB::Fasta', $genome;
	
print "Reading GFF...\n";
my $info = {};
open(IN, "$gff_file") or die ("Error reading $gff_file: $!\n");
while (<IN>){
  chomp;
  next if (/^\s*#/);
  my ($chr, $fucker, $feature_type, $f, $t,
      $score, $strand, $frame, $extra) = split(/\t/);
  $extra =~ s/\"//g;
  my @x = split(/\s*;\s*/, $extra);
  my $hextra = {};
  foreach my $x (@x){
    my ($k, $v) = split(/=/, $x);
    $hextra->{$k} = $v;
  }
  
  my $ID = $hextra->{'ID'};
  my $parent = $hextra->{'Parent'} || '';

  if ($feature_type eq 'gene'){
    $info->{$feature_type}{$ID}{'Name'} = $hextra->{'Name'};
    $info->{$feature_type}{$ID}{'Alias'} = $hextra->{'Alias'};
    $info->{$feature_type}{$ID}{'Note'} = $hextra->{'Note'};
    push @{$info->{'index'}}, $ID;
    $info->{$feature_type}{$ID}{'Son'} = $ID . '-RA';
  }
  if ($feature_type eq 'mRNA'){
    $info->{$feature_type}{$ID}{'Parent'} = $parent;
  }
  if ($feature_type eq 'exon'){
    $info->{$feature_type}{$ID}{'Chr'} = $chr;
    $info->{$feature_type}{$ID}{'From'} = $f;
    $info->{$feature_type}{$ID}{'To'} = $t;
    $info->{$feature_type}{$ID}{'Parent'} = $parent;
    push @{$info->{'mRNA'}{$parent}{'exon_son'}}, $ID;
  }
  if ($feature_type eq 'CDS'){
    $info->{$feature_type}{$ID}{'Chr'} = $chr;
    my $from = eval join '', $f, $strand, $frame;
 	my $to = eval join '', $t, $strand, $frame;
	if ($strand eq '-'){
      my $from = eval join '', $t, $strand, $frame;
 	  my $to = eval join '', $f, $strand, $frame;
	}
    $info->{$feature_type}{$ID}{'From'} = $from;
    $info->{$feature_type}{$ID}{'To'} = $to;
    $info->{$feature_type}{$ID}{'Parent'} = $parent;
	my $coors = $chr . ':' . $from . ',' . $to;
    push @{$info->{'mRNA'}{$parent}{'CDS_son'}}, $coors;
  }
}
close IN;

print "Rebuilding ORFs...\n";
system("mkdir seqs");
foreach my $i (@{$info->{'index'}}){
  my $fasta = '>';
  $fasta .= $info->{'gene'}{$i}{'Name'} . '; ' . $info->{'gene'}{$i}{'Alias'} . '; ' . $info->{'gene'}{$i}{'Note'} . "\n";
  my $seq;
  foreach my $j (@{$info->{'mRNA'}{$info->{'gene'}{$i}{'Son'}}{'CDS_son'}}){
    $seq .= $db{$j};
  }
  $seq =~ s/\n//g;
  $seq =~ s/(.{1,60})/$1\n/g;
  $fasta .= $seq;
  
  my $name = $info->{'gene'}{$i}{'Name'};
  
  my $outfile = 'NM_' . $info->{'gene'}{$i}{'Name'} . '.nfa';
  open (OUT, '>', $outfile) or die "Could not open $outfile, $!";
  print OUT $fasta;
  close OUT;
  
  system("sh gff2fasta_exe.sh $name");
}

print "Done!\n";

## Rebuild ORF
#if (exists($db{$chr})){
#  foreach my $e (@$tmp){
#    my $f = ($strand > 0) ? $e->{'from'} : $e->{'to'};
#    my $t = ($strand > 0) ? $e->{'to'} : $e->{'from'};
#    my $str = "$chr\:$f\,$t";
#    $seq .= $db{$str} if ($str);
#  }
#}

sub debug{
  my $hash = shift;
  my $i = 0;
  my $level = 0;
  my $MAX = 10;
  my $uid = "\&__$i\__";
  my $xml = '<root>'.$uid.'</root>';
  my $id_list = {$uid => $hash};
  while (scalar keys %$id_list){
    my $new_id_list = {};
    $level++;
    foreach my $id (keys %$id_list){
      my $temp_xml = '';
      my $href = $id_list->{$id};
      if (ref($href) eq 'ARRAY'){
        my $counter = 0;
        foreach my $val (@$href){
          $i++;
          $uid = "\&__$i\__";
          $new_id_list->{$uid} = $val;
          $temp_xml .= "\<c_$level\_$counter\>$uid\<\/c_$level\_$counter\>";
          $counter++;
          last if ($counter > $MAX);
        }
      }
      elsif (ref($href) eq 'HASH' || ref($href) eq __PACKAGE__){
        my $counter = 0;
        foreach my $key (keys %$href){
          $i++;
          $uid = "\&__$i\__";
          $new_id_list->{$uid} = $href->{$key};
          my $safe = '';
          if (substr($key,0,1) =~ /[^a-zA-Z]/){
            $safe = 'a';
          }
          $temp_xml .= "\<$safe$key\>$uid\<\/$safe$key\>";
          $counter++;
          last if ($counter > $MAX);
        }
      }
      else{
        $href =~ s/[<>]//g;
        $href = '<![CDATA['.$href.']]>';
        $temp_xml .= $href;
      }
      $temp_xml = '_empty_' unless ($temp_xml);
      die ("$id\t$temp_xml\n") unless ($xml =~ /$id/);
      $xml =~ s/$id/$temp_xml/;
    }
    $id_list = $new_id_list;
  }
  return $xml;
}