#!/usr/bin/perl -w
use strict;
use Cwd;
use LWP;
 
my $gene = shift || 'LMNA';
 
my $ua = LWP::UserAgent->new;
$ua->agent('Mozilla/5.0 (Windows NT 6.1; WOW64; rv:40.0) Gecko/20100101 Firefox/40.0');
$ua->parse_head(0);
 
my $datafile = 'ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current/homologene.data';
my $version = 'ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current/RELEASE_NUMBER';
 
my $localfile = getcwd() . '/' . 'homologdb.txt';
my $resultfile = getcwd() . '/' . $gene;
 
# Get ncbi version
my $v = remote_get($ua, $version);
chomp $v;
my $lv;
my @info;
if (-e $localfile){
  open(IN, $localfile);
  $lv = <IN>;
  chomp $lv;
  if ($v eq $lv){ # Use local file
    warn("Using local version $lv\n");
    @info = <IN>;
  }
  close IN;
}
if (!$lv || $lv ne $v){ # Download and use remote file
  warn("New remote version. Getting data...\n");
  my $response = remote_get($ua, $datafile);
  die("Could not get response for $datafile\n") unless($response);
  @info = split(/\n/, $response);
  my $tolocal = "$v\n" . $response;
  open(OUT, ">$localfile") or die("Cannot open $localfile: $!\n");
  print OUT $tolocal;
  close OUT;
}
 

warn("Loading data...\n");
my $homologs = {};
foreach my $l (@info){
  my ($gid, $taxid, $dumm, $symbol, $mygi, $protid) = split(/\s+/, $l);
  $homologs->{$symbol}{$taxid} = $mygi;
}
 
warn("Getting sequences...\n");
open (OUT, ">$resultfile") or die("Could not open $resultfile: $!\n");
foreach my $h (keys (%{$homologs->{$gene}})){
  warn("\tGetting $gene from $h...\n");
  my $url = 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?val=' .
            $homologs->{$gene}{$h} .
            '&db=protein&dopt=fasta&extrafeat=0&fmt_mask=0&maxplex=1&sendto=t&withmarkup=on&log$=seqview&maxdownloadsize=1000000';
 
  my $s = remote_get($ua, $url);
  my ($header, $seq) = $s =~ /^>([^\n]+)\n(.+)/sm;
  my @data = split(/\|/, $header);
  my $spinfo = pop @data;
  my $toshow = $spinfo;
  if ($spinfo =~ /\[([^\]]+)\]/){
    my @sp = split(/\s+/, $1);
    my $genus = shift @sp;
    my $sp = pop @sp;
    if ($genus && $sp){
      $toshow = substr($genus, 0, 1) . '_' . substr($sp, 0, 3);
    }
  }
  warn("\t$toshow\t->\t$spinfo\n");
  print OUT ">$toshow\n$seq\n";
}
close OUT;
 
 
 
sub remote_get{
  my $ua = shift;
  my $url = shift;
  my $result = '';
  my $tmp = $ua->get($url);
  if (exists($tmp->{'_content'})){
    $result = $tmp->{'_content'};
  }
  return $result;
}
 
sub hash_to_string {
  my $hash_ref = shift;
  my $tab      = shift;
  $tab = 1 unless ($tab);
  my $result = '';
  my @keys   = ();
  if (ref($hash_ref) eq 'HASH' ||
      ref($hash_ref) eq 'HTTP::Response' ||
      ref($hash_ref) eq 'HTTP::Headers' ||
      ref($hash_ref) eq 'HTTP::Request'
      ){
    @keys = sort { lc($a) cmp lc($b) } keys %$hash_ref;
  }
  elsif (ref($hash_ref) eq 'ARRAY'){
    @keys = sort { lc($a) cmp lc($b) } @$hash_ref;
  }
  foreach my $key ( @keys ) {
    my $key_val = '';
    if (ref($hash_ref) eq 'HASH' ||
        ref($hash_ref) eq 'HTTP::Response' ||
        ref($hash_ref) eq 'HTTP::Headers' ||
        ref($hash_ref) eq 'HTTP::Request'
    ){
      $key_val = $hash_ref->{$key};
    }
    if (ref($hash_ref) eq 'ARRAY'){
      $key_val = $key;
    }
    if ( ref( $key_val ) eq 'HASH' ) {
      my $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( \%{ $key_val }, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'HTTP::Headers' ) {
      my $text = HTTP::Headers->new;
      $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( $text, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'HTTP::Response' ) {
      my $text = HTTP::Response->new;
      $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( $text, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'HTTP::Request' ) {
      my $text = HTTP::Request->new;
      $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( $text, $tab + 1 );
    }
    elsif ( ref( $key_val ) eq 'ARRAY' ) {
      $result .= sp($tab)."$key => $key_val\n";
      $tab = $tab + 1;
      for (my $i = 0; $i < @{$key_val}; $i++){
        my $text = $key_val->[$i];
        my $t = sp($tab) . sp( length("$i => "), " " );
        $text =~ s/\n/\n$t/g;
        $result .=
            sp($tab)
          . "$i => $text\n";
        if (ref($key_val->[$i]) eq 'HASH'){
          $result .= hash_to_string( \%{ $key_val->[$i] }, $tab + 1 );
        }
        elsif (ref($key_val->[$i]) eq 'ARRAY'){
          $result .= hash_to_string( \@{ $key_val->[$i] }, $tab + 1 );
        }
      }
      $tab = $tab - 1;
    }
    else {
      my $text = $key_val;
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .= sp($tab) . "$key => $text\n";
    }
  }
  return $result;
}
 
sub sp {
  my $n   = shift;
  my $pad = shift;
  $pad = "\t" unless ($pad);
  return '' unless ($n);
  my $result;
  for ( my $i = 0 ; $i < $n ; $i++ ) {
    $result .= $pad;
  }
  return $result;
}
sub hash_to_xml {
 
  # Takes a hash and loads it into an xml-formatted string
  # @_ = (hash_ref, [tab], [previous_key])
  # hash_ref is a reference to a hash
  # tab is the indentation level
  # previous_key is the parent key to the current xml text
  # Returns an xml-formatted string with the information contained
  # in hash_ref.
  # This function is the reverse of get_xml(). The parameters
  # tab and previous_key are not used when calling the function.
  # The function itself uses them in recursive calls to itself.
 
  my $hash_ref = shift;
  my $tab      = shift;
  my $prev_key = shift;
  $tab = 1 unless ($tab);
  my $result = '';
  my @keys   = keys %$hash_ref;
  foreach my $key ( sort { lc($a) cmp lc($b) } @keys ) {
    if (ref($hash_ref->{$key}) eq 'HASH') {
      my $text = $hash_ref->{$key};
      if ($key =~ /^\d+$/){
        $result .=
          loadTag("$prev_key", hash_to_xml( \%{ $hash_ref->{$key} }, $tab, $prev_key ));
      }
      else{
        my @keys   = keys %{$hash_ref->{$key}};
        my $k = $keys[0];
        if (ref($hash_ref->{$key}{$k}) eq 'HASH' && $k =~ /^\d+$/){
          $result .= hash_to_xml( \%{ $hash_ref->{$key} }, $tab, $key )
        }
        else{
          $result .=
            loadTag( "$key", hash_to_xml( \%{ $hash_ref->{$key} }, $tab + 1, $key ));
        }
      }
    }
    elsif (ref($hash_ref->{$key}) eq 'HTTP::Headers'){
      my $text = HTTP::Headers->new;
      $text = $hash_ref->{$key};
      $result .=
            loadTag( "$key", hash_to_xml( $text, $tab + 1, $key ));
    }
    elsif (ref($hash_ref->{$key}) eq 'HTTP::Response'){
      my $text = HTTP::Response->new;
      $text = $hash_ref->{$key};
      $result .=
            loadTag( "$key", hash_to_xml( $text, $tab + 1, $key ));
    }
    elsif (ref($hash_ref->{$key}) eq 'HTTP::Request'){
      my $text = HTTP::Request->new;
      $text = $hash_ref->{$key};
      $result .=
            loadTag( "$key", hash_to_xml( $text, $tab + 1, $key ));
    }
    elsif ( ref( $hash_ref->{$key} ) eq 'ARRAY' ) {
      my $text = join (', ', @{$hash_ref->{$key}});
      $text =~ s/\&/\&amp;/g;
      $text =~ s/\"/\&quot;/g;
      $text =~ s/\'/\&apos;/g;
      $text =~ s/>/\&gt;/g;
      $text =~ s/</\&lt;/g;
      $result .= loadTag( "$key", $text );
    }
    else {
      my $text = $hash_ref->{$key};
      $text =~ s/\&/\&amp;/g;
      $text =~ s/\"/\&quot;/g;
      $text =~ s/\'/\&apos;/g;
      $text =~ s/>/\&gt;/g;
      $text =~ s/</\&lt;/g;
      $result .= loadTag( "$key", $text );
    }
  }
  return $result;
}
 
sub loadTag {
 
  # Loads information into a <tag>info</tag> structure
  # @_ = (tag_name, text)
  # Returns a scalar that contains the resulting structure
  # with newlines and tabs
  # The result can be used as <text> in another call to
  # create a nested structure
 
  my ( $tid, $ttext ) = @_;
  $ttext = '1' unless ($ttext);
  chomp $ttext;
  $ttext =~ s/(\n)/\n\t/gs;
  my $result = "<$tid>\n\t$ttext\n<\/$tid>\n";
  return $result;
}