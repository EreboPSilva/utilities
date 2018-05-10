#!/usr/bin/perl -w
use strict;

use Bio::SearchIO;

my $file  = shift;
my @temp = split /[\/\\]/, $0;
my $pname = pop @temp;
die("Use: perl $pname blast_file\n")
  unless ($file && -e $file);

my $blast = Bio::SearchIO->new('-format' => 'blast', '-file' => $file);


my $regions = {};
while(my $query = $blast->next_result()){
  my $seqname = $query->{'_queryacc'};
  my $hit = $query->next_hit(); # Only the best hit
  my $hitname = $hit->{'_hsps'}[0]{'-hit_name'};
  if ($hitname){
    my $hitdata = {
      'template' => $hit->{'_accession'},
      'expect'   => $hit->{'_significance'},
      'queryacc' => $seqname,
      'contig'   => $hitname,
      'score'    => $hit->{'_hsps'}[0]{'-score'},
      'from'     => $hit->{'_hsps'}[0]{'-hit_start'},
      'to'       => $hit->{'_hsps'}[0]{'-hit_end'},
    };
    push @{$regions->{$hitname}}, $hitdata;
  }
  #else{
  #  die(debug($hit, 50));
  #}
  #print "$seqname\t$hitname\n";
}

my $noverlap = {};
foreach my $region (keys %$regions){
  foreach my $hitdata (@{$regions->{$region}}){
    next if (exists($hitdata->{'used'}));
    my $result = [$hitdata];
    $hitdata->{'used'} = 1;
    foreach my $h2 (@{$regions->{$region}}){
      next if exists($h2->{'used'});
      if (overlap($hitdata, $h2)){
        $h2->{'used'} = 1;
        if ($h2->{'score'} > $hitdata->{'score'}){
          unshift @$result, $h2;
        }
        else{
          push @$result, $h2;
        }
      }
    }
    push @{$noverlap->{$region}}, $result;
  }
}

foreach my $region (sort{val($a) <=> val($b)} keys %$noverlap){
  foreach my $r (sort{$a->[0]{'from'} <=> $b->[0]{'from'}} @{$noverlap->{$region}}){
    my $acc = join("\t", map {$_->{'queryacc'}} @$r);
    my $l = abs($r->[0]{'to'} - $r->[0]{'from'});
    print "$region:$r->[0]{'from'}\-$r->[0]{'to'}\t$l\t$acc\n";
  }
}

sub val{
  my $v = shift;
  my $result = 0;
  if ($v =~ /(\d+)/){
    $result = $1;
  }
  return $result;
}

sub debug{
  my $hash = shift;
  my $i = 0;
  my $level = 0;
  my $MAX = shift || 10;
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
      elsif (ref($href) ){
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

sub overlap {

  # Checks whether two segments overlap
  # @_ = (hash_ref1, hash_ref2)
  # @_ = (number, hash_ref2)
  # @_ = (hash_ref1, number)
  # Returns 1 if there is overlapping, and 0 otherwise
  # If a hash is passed, it must have two keys called 'from'
  # and 'to' containing numerical values. These numbers define
  # a segment. The function succeeds if:
  #   a) Two segments from two passed hash_ref overlap (even
  # if at one point)
  #   b) A passed number is inside a passed segment (including
  # extremes)
  #   c) Two passed numbers are the same. This use would be
  # trivial and inelegant

  my $hash1 = shift;
  my $hash2 = shift;
  unless ( ref($hash1) eq 'HASH' ) {
    my %temp = ();
    $temp{'from'} = $hash1;
    $temp{'to'}   = $hash1;
    $hash1        = \%temp;
  }
  unless ( ref($hash2) eq 'HASH' ) {
    my %temp = ();
    $temp{'from'} = $hash2;
    $temp{'to'}   = $hash2;
    $hash2        = \%temp;
  }
  if ( exists( $hash1->{'from'} )
    && exists( $hash2->{'from'} )
    && exists( $hash1->{'to'} )
    && exists( $hash2->{'to'} ) )
  {
    my $f1     = $hash1->{'from'};
    my $f2     = $hash2->{'from'};
    my $t1     = $hash1->{'to'};
    my $t2     = $hash2->{'to'};
    my $result = 0;
    $result = 1 if ( ( $t2 - $t1 ) * ( $t2 - $f1 ) <= 0 );
    $result = 1 if ( ( $f2 - $t1 ) * ( $f2 - $f1 ) <= 0 );
    $result = 1 if ( ( $t1 - $t2 ) * ( $t1 - $f2 ) <= 0 );
    $result = 1 if ( ( $f1 - $t2 ) * ( $f1 - $f2 ) <= 0 );
    return $result;
  }
  else {
    return 0;
  }
}