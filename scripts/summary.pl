#!/usr/bin/perl -w

use strict;
use Pod::Usage;

# Very simple options manager

my %opts = ();
my %settings      = ();
my @done = ();

foreach my $arg (@ARGV){
  if ($arg =~ /\-(.+)=(.+)/){
    push @{$opts{'options'}}, $arg;
  }
  elsif (left($arg, 1) eq '-'){
    push @{$opts{'arguments'}}, $arg;
  }
  else{
    push @{$opts{'files'}}, $arg;
  }
}

if (exists($opts{'arguments'})){
  foreach my $arg (@{$opts{'arguments'}}){
    if ($arg eq '--help' || $arg eq '-h'){
      pod2usage(-verbose => 2);
    }
  }
}

my $settings_path = shift @{$opts{'files'}};
my $file          = shift @{$opts{'files'}};                   #File to edit
$settings_path = '' unless ($settings_path);

my $stag          = 1;
while ($stag) {
  $stag          = 0;
  %settings      = get_settings($settings_path, 'dbpath');
  while (!exists($settings{'dbpath'})){
    show_welcome();
    my $option = prompt();
    if ( $option eq '2' ) {
      my $file = load_project();
      %settings = get_settings($file, 'dbpath');
    }
    elsif ( $option eq '1' ) {
      %settings = prompt_for_settings();
    }
    elsif ( $option eq '3' ) {
      die("Bye!\n");
    }
    else {
      print "Please, enter the number of the option\n";
      print "or Ctrl-Z to exit\n";
    }
  }
  $settings_path = '';
  $stag          = 1 unless ( exists( $settings{'dbpath'} ) );
}
my @tbns   = get_files( $settings{'tbnpath'}, '.tbn$' );


# Check project results for inconsistencies

my $warnings = '';
my $results  = {};
my $hits     = {};
print "Reading...\n";
my $rawprot = '';
foreach my $tbn (@tbns){
  my $foldname = strip_file_name($tbn);
  my $folder   = "$settings{'basepath'}$foldname";
  my $tblastn  = "$settings{'tbnpath'}$tbn";
  my $sfile    = "$settings{'basepath'}$foldname\/$foldname\.seq";
  my $path     = "$settings{'basepath'}$foldname\/$foldname\.xml";
  my $info     = "$settings{'basepath'}$foldname\/info\.txt";
  my $wrn      = "$settings{'basepath'}$foldname\/warnings\.txt";
	my $protfile = "$settings{'basepath'}$foldname\/$foldname\.aa";
  if (-d $folder){
    my $inf = '';
    if (-e $info){
      $inf = get_file_text($info);
      if ($inf =~ /passed/i){
        $warnings .= "$foldname passed\n";
      }
    }
    unless ($inf){
      $warnings .= "No info.txt file in $foldname\n";
    }
    if ($inf && -e $path){
      if ($inf =~ /passed/i){
        $warnings .= "$foldname was passed, but xml file exists\n";
      }
      else{
        if (!get_file_text($wrn)){
          push @done, $foldname;
        }
        # Load results
        my $xml = get_file_text($path);
        my $gene = get_xml($xml);
				my $nseq = get_seq($sfile);
        numerate(\%{$gene->{'gene'}}, 'exon');
        my $temp = {};
        $temp = \%{$gene->{'gene'}{'exon'}};
        my $chr = $temp->{0}{'chromosome'};
        unless ($chr){
          $warnings .= "No chromosome field in $foldname\n";
        }
        else{
          foreach my $hit (keys %$temp){
            $temp->{$hit}{'name'} = $foldname;
            push @{$results->{$chr}}, $temp->{$hit};
          }
        }
				#Translate
				if ($nseq){
					my $prot = translate($nseq);
					$rawprot = $prot;
					$rawprot =~ s/[^\w\*]//g;
					open (OUT2, ">$protfile") or die("$protfile\: $!\n");
					print OUT2 ">$foldname\n";
					print OUT2 $prot;
					close OUT2;
					if ($rawprot =~ /\*./){
						$warnings .= "Early stop on $foldname\n";
					}
				}
        # Load tblastn hits
        if (-e $tblastn){
          my $ttext = get_file_text($tblastn);
          my %temp = loadBlast($ttext, $foldname);
          foreach my $chr (keys %temp){
            foreach my $hit (keys %{$temp{$chr}{'hits'}}){
              push @{$hits->{$chr}}, $temp{$chr}{'hits'}{$hit};
            }
          }
        }
        else{
          $warnings .= "$tblastn does not exist\n";
        }
      }
    }
  }
  else{
    $warnings .= "$folder does not exist. Please, run bsniffer on $foldname\n";
  }
}
$warnings .= "\n\nGenes without warnings\n";
$warnings .= join("\n", sort{$a cmp $b} @done);
print "Sorting...\n";
# Sort genes in each chromosome
foreach my $chr ( keys %$results ){
  my @temp = sort{$a->{'template'}{'from'}
                  <=>
                  $b->{'template'}{'to'}}
                  @{$results->{$chr}};
  $results->{$chr} = \@temp;
}

foreach my $chr (sort{ mynum($a) <=> mynum($b) || $a cmp $b } keys %$results){
  $warnings .= "\n\n\nChromosome $chr:\n\n";
  foreach my $exon ( @{$results->{$chr}} ){
    if (exists($exon->{'warnings'})){
      my $w = '';
      if (ref($exon->{'warnings'}) eq 'HASH'){
        $w = "Exon $exon->{'number'} at $exon->{'name'}\:".hjoin('. ', \%{$exon->{'warnings'}})."\n";
      }
      elsif ($exon->{'warnings'}){
        $w = "Exon $exon->{'number'} at $exon->{'name'}\: $exon->{'warnings'}\n";
      }
      $warnings .= $w;
    }
  }
}

open (OUT, ">summary.txt");
print OUT $warnings;
close OUT;
print "Done!";

sub get_seq{
	my $file = shift;
	my $result = '';
	if (-e $file){
		open (IN, $file);
		my $fst = <IN>;
		$result .= $fst unless ($fst =~ />/);
		while (<IN>){
			chomp;
			s/\W//g;
			$result .= $_;
		}
		close IN;
	}
	return $result;
}

sub mynum{
  my $text = shift;
  $text =~ s/\D//g;
  $text = '0' unless ($text);
  return $text;
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
  if ( !exists( $hash1->{'from'} ) ) {
    my %temp = ();
    $temp{'from'} = $hash1;
    $temp{'to'}   = $hash1;
    $hash1        = \%temp;
  }
  if ( !exists( $hash2->{'from'} ) ) {
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


sub get_file_name {

  #Takes a file path and returns the file name, with extension
  my $fpath = shift;
  my @path  = split( /\//, $fpath );
  my $l     = @path;
  my $rname = $path[ $l - 1 ];
  chomp $rname;
  return $rname;
}

sub strip_file_name {

  #Takes a file path and returns only the file name, without extension
  my $fpath = shift;
  $fpath = get_file_name($fpath);
  my @s = split( /\./, $fpath );
  return $s[0];
}

sub get_file_text{
  # Returns the whole file text in a single variable
	my $infile = shift;
	my $result = '';
	if (-e $infile){
  	my $backup = $/;
  	open (IN, $infile) or die ("Error while opening $infile: $!\n");
  	$/ = undef;
  	$result = <IN>;
  	close IN;
  	$/ = $backup;
	}
	return $result;
}


sub get_settings {

  # Gets settings from settings file and checks that certain keys are
  # present
  # @_ = (path, [mandatory_setting], [mandatory_setting], ...)
  # Returns a hash with each setting and each corresponding value

  my $path    = shift;
  my @tocheck = @_;
  my %result = ();
  $path = "projects\/$path"    if ( -e "projects\/$path" );
  $path = "projects\/$path.gt" if ( -e "projects\/$path.gt" );
  if ( $path && -e $path && !-d $path ) {
    open( IN, $path ) or die("Could not open $path\: $!");
    while (<IN>) {
      chomp;
      s/\s+=\s+/=/g;
      my ( $key, $val ) = split(/=|, |; /);
      if ( $key && $val){
        $result{$key} = $val;
      }
    }
    close IN;
    if (@tocheck){
      foreach my $checkme (@tocheck){
        unless ( exists( $result{$checkme} ) ) {
          print "Sorry, $path is not a correct project file\n";
          %result = ();
          return %result;
        }
      }
    }
  }
  else {
    %result = ();
    return %result;
  }
  return %result;
}

sub prompt_for_settings {
  my %result = ();
  print "***********                                       \n";
  print "New project                                       \n";
  print "***********                                       \n";
  print "                                                  \n";
  print "Please enter the settings to create a new project. Please\n";
  print "note that you need to create the folder in which you are \n";
  print "going to store your results before running this script.  \n";
  print "Always use '/' as a folder separator                     \n";
  print "To cancel, press Ctrl+C.                   \n";
  my $name = '';

  while ( !$name ) {
    $name = prompt('How do you want to name the new project?');
    unless ( -d 'projects' ) {
      mkdir 'projects';
    }
    if ( -e "projects\/$name.gt" ) {
      my $r = prompt("Project name already exists. Overwrite? \[y\/n\]");
      unless ( $r eq 'y' || $r eq 'yes' ) {
        $name = '';
        my @imposs = get_files( 'projects', '.gt' );
        print "Existing projects:\n";
        print join( "\n", @imposs ) . "\n";
      }
    }
  }
  $result{'basepath'} =
    ask( 'Where do you want your results stored?', 'folder' );
  $result{'dbpath'} =
    ask( 'Which file contains the template sequence?', 'file' );
  $result{'tbnpath'} = ask( 'Which path to the tblastn files?', 'folder' );
  $result{'ppath'} =
    ask( 'Which path to the known protein sequences?', 'folder' );
  $result{'aa_ext'} =
    prompt('Which is the extension of the protein sequence files?');
  $result{'npath'} =
    ask('Which path to the known nucleotide sequences? (optional)', 'folder');
  $result{'nt_ext'} =
    prompt(
    'Which is the extension of the nucleotide sequence files? (optional)');
  open( OUT, ">projects\/$name.gt" ) or die ("Could not open projects\/$name: $!\n");
  foreach my $key ( keys %result ) {
    print OUT "$key\=$result{$key}\n";
  }
  close OUT;
  return %result;
}
sub ask {

  # Prompts the user for a folder or a file.
  # @_ = (question_string, type_string)
  # question_string: prompt for the user
  # type_string: either 'file' or 'folder'
  # Returns a valid file or folder name.

  my $question    = shift;
  my $answer_type = shift;
  return '' unless ($answer_type);
  my $result = '';
  my $answer = '';
  while ( !$result ) {
    print "$answer is not a valid $answer_type\n" if ($answer);
    $answer = prompt($question);
    if ( $answer_type eq 'folder' && -d $answer ) {
      $answer .= '/' unless ( right( $answer, 1 ) eq '/' );
      $result = $answer;
    }
    if ( $answer_type eq 'file' && -e $answer && !-d $answer ) {
      $result = $answer;
    }
  }
  return $result;
}

sub load_project {
  my $tag    = 1;
  my $result = '';
  while ($tag) {
    $tag = 0;
    display_saved_projects();
    my $fp = prompt("Which project do you want to open?");
    $fp = "projects\/$fp";
    if ( $fp && ( -e $fp || -e "$fp.gt" ) ) {
      if ( -e $fp ) {
        $result = $fp;
      }
      elsif ( -e "$fp.gt" ) {
        $result = "$fp.gt";
      }
    }
    else {
      print "$fp does not exist!\n\n";
      $tag = 1;
    }
  }
  return $result;
}

sub show_welcome {
  print "\n\n****************************************************\n";
  print
"Welcome to GeneTuner, a tool to manually curate gene\npredictions by similarity\n\n";
  print "Good luck!\n\n\n";
  print "Please, choose an option:\n\n";
  print "\t1.- Create a new project\n\n";
  print "\t2.- Load a project\n\n";
  print "\t3.- Quit\n\n";
}

sub display_saved_projects {
  print "Here is a list of saved projects:\n";
  my @files = get_files( 'projects', ".gt\$" );
  print join( "\n", @files ) . "\n";
}

sub prompt {
  my $disp = shift;
  print "$disp\n" if ($disp);
  my $input = <STDIN>;
  chomp $input;
  return $input;
}

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

sub get_files {

  # Gets the files in a directory
  # @_ = ([folder_path], [filter])
  # Returns an array with the names of all the files in folder_path
  # that comply with filter sorted by alphabetical order
  # If present, filter must be a regex
  # No folder_path means current folder. Cannot be omited if
  # a filter is passed
  # No filter means all files

  my $dpath  = $_[0] ? $_[0] : ".";
  my $filter = $_[1] ? $_[1] : "";
  my @result = ();
  my @unsort = ();
  $filter =~ s/\\/\\\\/gs;    # Interpolate
  opendir DIR, $dpath or die "unable to open $dpath\n";
  while ( my $name = readdir(DIR) ) {
    if ( ( $name =~ /$filter/ ) ) {
      push @unsort, $name;
    }
  }
  closedir DIR;
  @result = sort{ $a cmp $b } @unsort;
  return @result;
}

sub get_xml {

  # Takes an xml-formatted string and loads it into a hash
  # @_ = (xml_text)
  # Returns a ref to a hash following the structure of the xml text.
  # Tags are turned into keys, whose value is what is inside the tag.
  # If a key is repeated, a subhash is created for each instance. The
  # key for each subhash is its correlative number.
  # Limitations: if a tag contains text and other tags, the text is
  # ignored. Modifiers inside the tag would cause failures.

  my $xml  = shift;
  my $hash = shift;
  $hash = {} unless ($hash);
  while ( $xml =~ /<(.+?)>(.+?)<\/\1>/gs ) {
    my $tag     = $1;
    my $content = $2;
    $content =~ s/[\s\n]+$//;
    $content =~ s/^[\s\n]+//;
    if ( $content =~ /<(.+?)>(.+?)<\/\1>/s ) {
      if ( exists( $hash->{$tag} ) && !exists( $hash->{$tag}{0} ) ) {
        my $vzero = $hash->{$tag};
        delete $hash->{$tag};
        $hash->{$tag}{0} = $vzero;
        $hash->{$tag}{1} = get_xml( $content, \%{ $hash->{$tag}{1} } );
      }
      elsif ( exists( $hash->{$tag}{0} ) ) {
        my $n = scalar keys %{ $hash->{$tag} };
        $hash->{$tag}{ $n } =
          get_xml( $content, \%{ $hash->{$tag}{ $n } } );
      }
      else {
        $hash->{$tag} = get_xml( $content, \%{ $hash->{$tag} } );
      }
    }
    else {
      if ( exists( $hash->{$tag} ) && ref($hash->{$tag}) ne 'HASH' ) {
        my $vzero = $hash->{$tag};
        delete $hash->{$tag};
        $hash->{$tag}{0} = $vzero;
        $hash->{$tag}{1} = $content;
      }
      elsif (  ref($hash->{$tag}) eq 'HASH') {
        my $n = scalar keys %{$hash->{$tag}};
        $hash->{$tag}{ $n } = $content;
      }
      else{
        $hash->{$tag} = $content;
      }
    }
  }
  return $hash;
}

sub numerate{
  # Transforms a subhash into a cumulative subhash. Useful to correct the behaviour
  # of get_xml.
  my $temp = shift;
  my $key  = shift;
  if (exists($temp->{$key})){
    my $temp2 = {};
    unless (ref($temp->{$key}) eq 'HASH' && exists($temp->{$key}{0})){
      $temp2->{0} = $temp->{$key};
    }
    else{
      $temp2 = $temp->{$key};
    }
    $temp->{$key} = $temp2;
  }
}

sub loadBlast {

# Warning: This version of loadBlast sepparates exons inside a hit
# @_ = (blast_text)
# returns a hash containing info about each blast hit sorted from lower to higher
# by the 'from' field:
#
#	contig_number =>
#		contig_name  => (contig_name)
#		hits =>
#			hit_number =>--------------------|
#				expect => (expect)             |
#				from   => (template_from)      |
#				ids    => (ids, 0-1)           |
#				query  =>                      |
#					from   => (query_from)       |
#					seq    => (query_seq)        |-> Returned by get_hit_data()
#					to     => (query_to)         |
#				score  => 49                   |
#				template =>                    |
#					from   => (template_from)    |
#					seq    => (template_seq)     |
#					to     => (template_to)      |
#         strand => (strand)           |
#				to => (template_to)------------|
#
# Example:
#
# >chr10_random
#           Length = 39000292
#
#  Score = 49.3 bits (116), Expect = 6e-04
#  Identities = 30/77 (38%), Positives = 44/77 (57%), Gaps = 6/77 (7%)
#  Frame = -2
#
# Query: 400      VPGYVSFRHVEEGSWYIQSLCNHLKKLVPRHEDILSILTAVNDDVSRRVDKQG------T 453
#                 V GY S+R    GSW++Q+LC+ L++   +  +I+ ILT VND V+R  + Q
# Sbjct: 19893264 VEGYYSWRSPGRGSWFVQALCSILEE-HGKDLEIMQILTRVNDRVARHFESQSDDPRFHE 19893088
#
# Query: 454      KKQMPQPAFTLRKKLVF 470
#                 KKQ+P     L K+L F
# Sbjct: 19893087 KKQIPCVVSMLTKELYF 19893037
#
#
#
#  Score = 39.7 bits (91), Expect = 0.49
#  Identities = 24/77 (31%), Positives = 38/77 (49%), Gaps = 11/77 (14%)
#  Frame = +2
#
# Query: 357      ACQGEEIQPSVSIEADALNPEQAPTSLQDSIPAEADFLLGLATVPGYVS---------FR 407
#                 AC+G E+   +  ++  +N   A    +  IP EADFL   +TVPG  S         FR
# Sbjct: 22007009 ACRGTELDDGIQADSGPINDTDANPRYK--IPVEADFLFAYSTVPGITSIVCQAYFSPFR 22007182
#
# Query: 408      H--VEEGSWYIQSLCNH 422
#                 H   + G+ ++ S C++
# Sbjct: 22007183 HQNTDLGTVFVSSCCSN 22007233
#
# The result would be:
#	0 => (contig_number)
#		contig_name => chr10_random
#		hits =>
#			0 => (hit_number)
#				expect => 6e-04
#				from => 19893264
#				ids => 0.38961038961039
#				query =>
#					from => 400
#					seq => VPGYVSFRHVEEGSWYIQSLCNHLKKLVPRHEDILSILTAVNDDVSRRVDKQG------TKKQMPQPAFTLRKKLVF
#					to => 470
#				score => 49
#				template =>
#					from => 19893264
#					seq => VEGYYSWRSPGRGSWFVQALCSILEE-HGKDLEIMQILTRVNDRVARHFESQSDDPRFHEKKQIPCVVSMLTKELYF
#					to => 19893037
#				to => 19893037
#			1 =>
#				expect => 0
#				from => 22007009
#				ids => 0.311688311688312
#				query =>
#					from => 357
#					seq => ACQGEEIQPSVSIEADALNPEQAPTSLQDSIPAEADFLLGLATVPGYVS---------FRH--VEEGSWYIQSLCNH
#					to => 422
#				score => 39
#				template =>
#					from => 22007009
#					seq => ACRGTELDDGIQADSGPINDTDANPRYK--IPVEADFLFAYSTVPGITSIVCQAYFSPFRHQNTDLGTVFVSSCCSN
#					to => 22007233
#				to => 22007233
#

  my $btext   = shift;
  my $tag     = shift;
  my $header  = '';
  my $foot    = '';
  my %result  = ();
  my @contigs = isplit( '>', $btext );
  $header = shift @contigs unless ( substr( $contigs[0], 0, 1 ) eq '>' );
  my $ncontigs = @contigs;
  for ( my $i = 0 ; $i < $ncontigs ; $i++ ) {
    my $contig_text = $contigs[$i];
    $contig_text =~ />(.+)\n/;
    my $cname         = $1;
    my %temp_hit_data = ();
    my %hit_data      = ();
    my @hits          = isplit( 'Score', $contig_text );
    foreach my $hit (@hits) {
      if ( lc( substr( $hit, 0, length('Score') ) ) eq 'score' ) {
        my %temp_hit = get_hit_data($hit, $tag);
        hpush( \%{ $hit_data{'hits'} }, \%temp_hit );
      }
    }
    $result{$cname} = \%hit_data;
  }
  return %result;
}

sub get_hit_data {

  # Describes a blast hit (see loadBlast)
  # @_ = (hit_text)
  # Returns a hash containing data from the hit_text
  # See loadBlast for an example

  my $text        = shift;
  my $tag         = shift;    #Text tag inserted in every hit
  my %result      = ();
  my $getn        = '.*?([-\\+]?\d?e?-?[\\d]+)';
  my $get_score   = ".*?Score$getn";
  my $get_expect  = ".*?Expect$getn";
  my $get_ids     = ".*?Ident$getn$getn$getn";
  my $get_seq     = '\s([^0-9]+)\s';
  my $align_query = "Query$getn$get_seq$getn.*?Sbjct$getn$get_seq$getn";
  my $total_query = $align_query;
  $total_query =~ s/[\(\)]//g;
  $text        =~ /$get_score$get_expect$get_ids/si;
  my $score  = $1;
  my $expect = $2;
  if ($expect && left($expect, 1) eq 'e'){
    $expect = "1$expect";
  }
  my $ids    = $5;
  $ids = $3 / $4 if ( $3 && $4 );

  my @lines =
    isplit( "Query$getn$get_seq$getn.*?Sbjct$getn$get_seq$getn", $text );
  my $nl     = @lines;
  my $totext = '';
  my $qtext  = '';
  my $ttext  = '';
  my $qfrom  = 0;
  my $tfrom  = 0;
  my $qto    = 0;
  my $tto    = 0;
  if ($nl) {

    for ( my $i = 1 ; $i < $nl ; $i++ ) {
      my $line = $lines[$i];
      $line =~ /$align_query/si;
      $qfrom = $1 if ( $i == 1 );
      $qto   = $3 if ( $i == $nl - 1 );
      $qtext .= $2;
      $tfrom = $4 if ( $i == 1 );
      $tto   = $6 if ( $i == $nl - 1 );
      $ttext .= $5;
      $line =~ /($total_query)/s;
      $totext .= "$1\n\n";
    }
    $qtext =~ s/[^-*\w]//g;
    $ttext =~ s/[^-*\w]//g;
    $result{'text'}   = $totext;
    $result{'score'}  = $score;
    $result{'expect'} = $expect;
    $result{'ids'}    = $ids;

    if ( $tfrom < $tto ) {
      $result{'strand'} = 1;
    }
    else {
      $result{'strand'} = -1;
    }

    $result{'from'}   = $tfrom;# - 1;
    $result{'to'}     = $tto;# - 1;
    $result{'tag'}    = $tag;
#    $result{'template'}{'from'} = $tfrom;# - $result{'strand'};
#    $result{'template'}{'to'}   = $tto;# - $result{'strand'};
#    $result{'template'}{'seq'}  = $ttext;
#    $result{'query'}{'from'}    = $qfrom;
#    $result{'query'}{'to'}      = $qto;
#    $result{'query'}{'seq'}     = $qtext;

    #$result{'qtext'} = delinsert($qtext, $ttext);
  }

#%result = () if ($result{'to'} < $result{'from'});  # Empty if the frame is negative
  return %result;
}

sub hpush {

  # Pushes a value into a hash with a correlative number as key
  # @_ = (hash_ref, value)
  # Returns the index of the entered value
  # Used in the context of cumulative hashes. The value can in turn
  # be a hash_ref to construct HoHs.

  my $hash  = shift;
  my $value = shift;
  my $n     = scalar keys %$hash;
  $hash->{$n} = $value;
  return $n;
}

sub hadd {

  # Adds a key->value pair inside the last key pushed to the target hash
  # @_ = (hash_ref, key, value)
  # No return value
  # Used in the context of cumulative hashes.

  my $hash  = shift;
  my $key   = shift;
  my $value = shift;
  my $n     = scalar keys %$hash;
  $hash->{ $n - 1 }{$key} = $value;
}

sub isplit {

  # Split that includes the pattern in the result
  #
  # @_ = (pattern, text);
  #
  # Returns an array containing <text> split by <pattern>,
  # without deleting each occurence of <pattern> in <text>, as
  # split does
  #
  # Example: isplit('\d', "7string3split5by6numbers8") =
  #             (7string, 3split, 5by, 6numbers, 8)
  # , whereas
  #           split(/\d/, "7string3split5by6numbers8") =
  #             (, string, split, by, numbers);

  my $re   = shift;
  my $text = shift;
  my @hits = ();
  my $l    = length($text);
  my $ppos = pos($text) = 0;
  while ( $text =~ /($re)/gs ) {
    my $cpos = pos($text) - length($1);
    if ( $cpos > $ppos ) {
      my $stext = substr( $text, $ppos, $cpos - $ppos );
      push @hits, $stext;
      $ppos = $cpos;
    }
  }
  push @hits, substr( $text, $ppos, $l - $ppos );
  return @hits;
}

sub hjoin{
  my $sep  = shift;
  my $hash = shift;
  my $temp = {};
  my $result = '';
  foreach my $key (%$hash){
    if ($key =~ /^\d+$/){
      $temp->{$key} = $hash->{$key};
    }
  }
  foreach my $key (sort{$a <=> $b} keys %$temp){
    if ($result){
      $result .= "$sep$temp->{$key}";
    }
    else{
      $result = $temp->{$key};
    }
  }
  return $result;
}

sub translate {
  my $seq   = shift;
  my %cdict = (
    'TTT' => 'F', 'TTC' => 'F', 'TTA' => 'L', 'TTG' => 'L',
    'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
    'TAT' => 'Y', 'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*',
    'TGT' => 'C', 'TGC' => 'C', 'TGA' => '*', 'TGG' => 'W',
    'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
    'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
    'CAT' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
    'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
    'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'ATG' => 'M',
    'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
    'AAT' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
    'AGT' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
    'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
    'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
    'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
    'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',
    'TCN' => 'S', 'CGN' => 'R', 'CCN' => 'P', 'CTN' => 'L',
    'ACN' => 'T', 'GGN' => 'G', 'GCN' => 'A', 'GTN' => 'V',
		
		'TTY' => 'F', 'TTR' => 'L',
		'TCY' => 'S', 'TAR' => '*',
		'CAY' => 'H', 'CAR' => 'Q',
		'TGY' => 'C',
		'ATY' => 'I', 'ATM' => 'I', 'ATW' => 'I', 'ATH' => 'I',
		'AAY' => 'N', 'AAR' => 'K',
		'AGY' => 'S', 'AGR' => 'R',
		'GAY' => 'D', 'GAR' => 'E',
		'MGA' => 'R', 'MGG' => 'R', 'MGR' => 'R',
		'YTA' => 'L', 'YTG' => 'L', 'YTR' => 'L'		
  );
	foreach my $d (qw/R Y S W K M B D H V/){
		$cdict{'TC'.$d} = 'S';
		$cdict{'CG'.$d} = 'R';
		$cdict{'CC'.$d} = 'P';
		$cdict{'CT'.$d} = 'L';
		$cdict{'AC'.$d} = 'T';
		$cdict{'GG'.$d} = 'G';
		$cdict{'GC'.$d} = 'A';
		$cdict{'GT'.$d} = 'V';
	}
  my $result = '';
  my @codons = unpack( '(A3)*', uc($seq) );
  foreach my $codon (@codons) {
    if ( exists( $cdict{$codon} ) ) {
      $result .= "$cdict{$codon}";
    }
    else {
      $result .= "X";
    }
  }
	my @t = $result =~ /(.{0,60})/g;
  $result = join("\n", @t);
  return $result;
}

################# Debug ###################
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
