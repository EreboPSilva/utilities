#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use Pod::Usage;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Storable;
use vBIO::BATI;
#require Term::Screen;
#use Win32::Console;
use strict;

my $terminal;

eval 'require Term::Screen';
if ($@){
  eval 'use Win32::Console';
  if ($@){
    die("You need either Term::Screen in *nix systems or Win32::Console in Windows\n");
  }
  else{
    $terminal = GT_Win32_Screen->new();
  }
}
else{
  $terminal = GT_Term_Screen->new();
}

my $bftr = vBIO::BATI->new('GeneTuner');

####Global variables####
my $FOLDER_SEP = "\/";      #Folder separation symbol
my $SPLICE_OVERHANG = 2;    #Number of positions checked for each splice site
my $CONTIG_OVERHANG =
  50000;    #Number of base pairs read from the chromosome at either
           #side of the gene
my $MAJOR_STEP    = 1000;                    #Length of 'fast scroller' step
#my $scr           = '';                      #Term::Screen object is stored here
#my $scr_in        = '';                      #Win32::Console object is stored here
#my $scr_out       = '';                      #Win32::Console object is stored here
my $MAP_WIDTH     = 60;                      #Line width in HTML map
my $MAP_OVERHANG  = 100;                     #Overhang in HTML map
my $INTRON_COLLAPSE = 5;                     #Intron can collapse in HTML map
                                             #if it spans more than $I_C lines
my %settings      = ();
my %db            = ();                      #Stores Bio::DB::Fasta object
my $CURRENT_INFO  = '';
my $ctrlz  = chr(26);
my $ctrlf  = chr(6);
my $ctrlb  = chr(2);
my $ctrlg  = chr(7);
my $enter  = chr(13);
########################

# Very simple options manager

my %opts = ();

foreach my $arg (@ARGV){
  if ($arg =~ /\-(.+)=(.+)/){
    push @{$opts{'options'}}, $arg;
  }
  elsif ($bftr->left($arg, 1) eq '-'){
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
$file          = '' unless ($file);

my $stag          = 1;
while ($stag) {
  $stag          = 0;
  %settings      = $bftr->get_settings($settings_path, 'dbpath');
  while (!exists($settings{'dbpath'})){
    $bftr->show_welcome();
    my $option = $bftr->prompt();
    if ( $option eq '2' ) {
      my $file = $bftr->load_project();
      %settings = $bftr->get_settings($file, 'dbpath');
    }
    elsif ( $option eq '1' ) {
      %settings = $bftr->prompt_for_settings();
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
my @tbns   = $bftr->get_files( $settings{'tbnpath'}, "\.tbn\$" );
my $nfiles = @tbns;
my $fopen   = '';
my $counter = 0;
if ($file){
  $CURRENT_INFO = $bftr->path($settings{'basepath'}, $file, "info.txt");
}
else{
  while ( !$fopen && $counter < $nfiles ) {
    my $f = $tbns[$counter];
    my $tname = $bftr->strip_file_name($f);
    $CURRENT_INFO = $bftr->path($settings{'basepath'}, $tname, "info.txt");
    unless ( -e $CURRENT_INFO ){
      $fopen = 'y';
      $file = $f;
    }
    else{
      my $t = get_file_text( $CURRENT_INFO);
      unless ($t){
        $fopen = 'y';
        $file = $f;
      }
    }
    $counter++;
  }
}
my %vs = ();

die("no new 'tbn' files found!\n") unless ($file);
$vs{'file'} = $bftr->strip_file_name($file);
$vs{'up_overhang'}   = 1;
$vs{'down_overhang'} = 1;
#Lock gene for editing
open (OUT, ">$CURRENT_INFO") or
      die ("Problem writing $CURRENT_INFO: $!\n");
print OUT "Locked";
close OUT;
initVS( \%vs );
$terminal->initScr();
$terminal->showVS( \%vs );
my $ch = '';
while ( $ch ne 'q' ) {
  $ch = $terminal->get_char();
  if ($ch) {
    $ch = lc($ch);
    my $refresh = 1;
    if (exists($vs{'lock'}) && $vs{'lock'}){
      unless ($vs{'lock'} =~ /\s$ch\s/){
        $ch = '';
        $terminal->message(\%vs, 'Command blocked at this time');
        $terminal->showVS(\%vs);
      }
    }
    if ( $ch eq 'q' ) {
      ########       Quit      ########
      my $answer = lc( $terminal->tprompt( \%vs, "Save changes?", 'yn' ) );
      if ( $answer eq 'y' ) {
        $terminal->endScr();
        print "\nSetting splice frames. This may take a few minutes...\n";
        set_splice_frames( \%vs );
        print "Saving results...\n";
        savegene( \%vs );
      }
      else{
        open (OUT, ">$CURRENT_INFO") or
              die ("Problem writing to disk: $!\n");
        print OUT '';
        close OUT;
      }
    }
    else {
      ######## Process command ########
      if ( $ch eq 'kr' ) {
        $vs{'cpos'} += $vs{'step'};
      }
      elsif ( $ch eq 'kl' ) {
        $vs{'cpos'} -= $vs{'step'};
      }
      elsif ( $ch eq 'p' ) {
        $vs{'cpos'} += $MAJOR_STEP;
      }
      elsif ( $ch eq 'o' ) {
        $vs{'cpos'} -= $MAJOR_STEP;
      }
      elsif ( $ch eq '<' or $ch eq '>' ){
        my $cpos = realpos(\%vs, $vs{'cpos'});
        my $where = 'upstream';
        if ( $ch eq '<'){
          $vs{'up_overhang'} += 1;
        }
        if ( $ch eq '>'){
          $vs{'down_overhang'} += 1;
          $where = 'downstream';
        }
        $terminal->endScr();
        initVS( \%vs );
        $terminal->initScr();
        $vs{'cpos'} = localpos(\%vs, $cpos);
        $terminal->message( \%vs, "Added $CONTIG_OVERHANG bases $where" );
      }
      elsif ( $ch eq 's' ) {
        my $prev = -1;
        my $next = -1;
        my $p = prevfeat( \%{ $vs{'features'}{'exon'} }, $vs{'cpos'} );
        my $n = nextfeat( \%{ $vs{'features'}{'exon'} }, $vs{'cpos'} );
        $prev = prevfeat( \%{ $vs{'features'}{'tblastn'}{0}{'hits'} }, $vs{'features'}{'exon'}{$p}{'to'} )
          if ($p > -1);
        $next = nextfeat( \%{ $vs{'features'}{'tblastn'}{0}{'hits'} }, $vs{'features'}{'exon'}{$n}{'from'} )
          if ($n > -1);
        my $qfrom = 1;
        my $qto   = $vs{'qlength'};
        my $tfrom = 1;
        my $tto   = $vs{'cpos'};
        $qfrom = $vs{'features'}{'tblastn'}{0}{'hits'}{$prev}{'query'}{'to'}
            if ($prev > -1);
        $qto   = $vs{'features'}{'tblastn'}{0}{'hits'}{$next}{'query'}{'from'}
            if ($next > -1);
        $tfrom = $vs{'features'}{'tblastn'}{0}{'hits'}{$prev}{'to'}
            if ($prev > -1);
        $tto   = $vs{'features'}{'tblastn'}{0}{'hits'}{$next}{'from'}
            if ($next > -1);
        $vs{'settings'}{'local'}{'contig'}{'tblastn'} = {
          'query' => {
            'from' => $qfrom,
            'to'   => $qto
          },
          'template' => {
            'from' => $tfrom,
            'to'   => $tto
          }
        };
        # Optimization for short sequences
        $vs{'settings'}{'local'}{'tblastn'} = {
          'e' => $vs{'settings'}{'tblastn'}{'e'},
          'M' => 'PAM30',
          'G' => 9
        };
        localize(\%vs, 'tblastn', 1);
      }
      elsif ( $ch eq 'y' ){
        if (!exists($vs{'settings'}{'local'}{'tblastn'}{'e'}) ||
            $vs{'settings'}{'local'}{'tblastn'}{'e'} == 0){
          $vs{'settings'}{'local'}{'tblastn'}{'e'} = $vs{'settings'}{'tblastn'}{'e'};
        }
        $vs{'settings'}{'local'}{'tblastn'}{'e'} *= 10;
        localize(\%vs, 'tblastn');
      }
      elsif ( $ch eq 'h' ){
        if (!exists($vs{'settings'}{'local'}{'tblastn'}{'e'}) ||
            $vs{'settings'}{'local'}{'tblastn'}{'e'} == 0){
          $vs{'settings'}{'local'}{'tblastn'}{'e'} = $vs{'settings'}{'tblastn'}{'e'};
        }
        $vs{'settings'}{'local'}{'tblastn'}{'e'} /= 10;
        localize(\%vs, 'tblastn');
      }
      elsif ( $ch eq '(' ){
        $vs{'settings'}{'local'}{'contig'}{'tblastn'}{'template'}{'from'} = $vs{'cpos'};
        my $qfrom = $vs{'settings'}{'local'}{'contig'}{'tblastn'}{'query'}{'from'};
        my $prompt = "Starting position in the query? \($qfrom\)";
        $qfrom = $terminal->prompt_for_expression(\%vs, $prompt);
        if ($qfrom > 0 && $qfrom < $vs{'qlength'}){
          $vs{'settings'}{'local'}{'contig'}{'tblastn'}{'query'}{'from'} = $qfrom;
        }
      }
      elsif ( $ch eq ')' ){
        $vs{'settings'}{'local'}{'contig'}{'tblastn'}{'template'}{'to'} = $vs{'cpos'};
        my $qto = $vs{'settings'}{'local'}{'contig'}{'tblastn'}{'query'}{'to'};
        my $prompt = "Ending position in the query? \($qto\)";
        $qto = $terminal->prompt_for_expression(\%vs, $prompt);
        if ($qto > 0 && $qto < $vs{'qlength'}){
          $vs{'settings'}{'local'}{'contig'}{'tblastn'}{'query'}{'to'} = $qto;
        }
      }
      elsif ( $ch eq 'kd' ) {
        $vs{'cpos'} = getnext( \%{ $vs{'features'}{'exon'} }, $vs{'cpos'} );
      }
      elsif ( $ch eq 'ku' ) {
        $vs{'cpos'} = getprev( \%{ $vs{'features'}{'exon'} }, $vs{'cpos'} );
      }
      elsif ( $ch eq '.' ) {
        $vs{'cpos'} =
          getnext( \%{ $vs{'features'}{'tblastn'}{0}{'hits'} }, $vs{'cpos'} );
      }
      elsif ( $ch eq ',' ) {
        $vs{'cpos'} =
          getprev( \%{ $vs{'features'}{'tblastn'}{0}{'hits'} }, $vs{'cpos'} );
      }
      elsif ( $ch eq 'l' ) {
        $vs{'cpos'} =
          getnext( \%{ $vs{'features'}{'blastn'}{0}{'hits'} }, $vs{'cpos'} );
      }
      elsif ( $ch eq 'k1' ) {
        debug(hash_to_string(\%vs));
        $terminal->message(\%vs, "Added internal dump to \'debug.txt\'");
      }
      elsif ( $ch eq 'w' ){
        my $c = $vs{'current_exon'};
        if ( $c == -1){
          $terminal->message( \%vs, 'No exon selected for warnings');
        }
        elsif ( $c > -1){
          my $w = $terminal->prompt_for_expression(\%vs, 'Warning?');
          chomp $w;
          if ($w){
            if ($bftr->left($w, 1) ne '_'){
              $vs{'features'}{'exon'}{$c}{'warnings'} = {};
            }
            else{
              $w = $bftr->dleft($w, 1);
            }
            $bftr->hpush(\%{$vs{'features'}{'exon'}{$c}{'warnings'}}, $w);
          }
          else{
            delete($vs{'features'}{'exon'}{$c}{'warnings'});
            $terminal->message( \%vs, "No warnings");
          }
        }
      }
      elsif ( $ch eq 'k' ) {
        $vs{'cpos'} =
          getprev( \%{ $vs{'features'}{'blastn'}{0}{'hits'} }, $vs{'cpos'} );
      }
      elsif ( $ch eq 'z' ) {
        my $n = nextfeat( \%{ $vs{'features'}{'exon'} }, $vs{'cpos'} );
        if ( $n > -1 ) {
          setundo( \%vs );
          my $dist = $vs{'features'}{'exon'}{$n}{'to'} - $vs{'cpos'};
          $dist = 3 * int( $dist / 3 ) + 2;
          $vs{'features'}{'exon'}{$n}{'from'} =
            $vs{'features'}{'exon'}{$n}{'to'} - $dist;
        }
      }
      elsif ( $ch eq 'x' ) {
        my $n = prevfeat( \%{ $vs{'features'}{'exon'} }, $vs{'cpos'} );
        if ( $n > -1 ) {
          setundo( \%vs );
          my $dist = -$vs{'features'}{'exon'}{$n}{'from'} + $vs{'cpos'};
          $dist = 3 * int( $dist / 3 ) + 2;
          $vs{'features'}{'exon'}{$n}{'to'} =
            $vs{'features'}{'exon'}{$n}{'from'} + $dist;
        }
      }
      elsif ( $ch eq 'f' ) {
        $vs{'settings'}{'blastn'}{'e'} /= 10;
        refreshVS( \%vs );
        $terminal->message( \%vs, "Blastn Expect = $vs{'settings'}{'blastn'}{'e'}");
      }
      elsif ( $ch eq 'r' ) {
        $vs{'settings'}{'blastn'}{'e'} *= 10;
        refreshVS( \%vs );
        $terminal->message( \%vs, "Blastn Expect = $vs{'settings'}{'blastn'}{'e'}");
      }
      elsif ( $ch eq 'd' ) {
        $vs{'settings'}{'tblastn'}{'e'} /= 10;
        refreshVS( \%vs );
        $terminal->message( \%vs, "TBlastn Expect = $vs{'settings'}{'tblastn'}{'e'}");
      }
      elsif ( $ch eq 'e' ) {
        $vs{'settings'}{'tblastn'}{'e'} *= 10;
        refreshVS( \%vs );
        $terminal->message( \%vs, "TBlastn Expect = $vs{'settings'}{'tblastn'}{'e'}");
      }
      elsif ( $ch eq 'm' ) {
        $vs{'cpos'} = getprevM( \%vs );
      }
      elsif ( $ch eq 'a' ) {
        setundo( \%vs );
        delexon( \%vs );
      }
      elsif ( $ch eq ' ' ) {
        if ( $vs{'step'} == 3 ) {
          $vs{'step'} = 1;
          delete $vs{'lock'};
        }
        else {
          setundo( \%vs );
          insertexon( \%vs );
          $vs{'lock'} = " kl kr ku kd \. \, l k e d   o p $ctrlg ";
        }
      }
      elsif ( $ch eq $ctrlz ) {    #Ctrl-Z
        if ( exists( $vs{'undo'} ) ) {
          %{ $vs{'features'}{'exon'} } =
            %{ Storable::dclone( \%{ $vs{'undo'} } ) };
        }
      }
      elsif ( $ch eq $ctrlf ) {
        my $regex = get_regex(\%vs);
        if ($regex){
          $terminal->message( \%vs, 'Press F3 to search downstream or F2 to search upstream' );
        }
      }
      elsif ( $ch eq $ctrlb ) {
        my $regex = get_regex(\%vs);
        if ($regex){
          $vs{'search_pattern'} = backtranslate($regex);
          $terminal->message( \%vs, 'Press F3 to search downstream or F2 to search upstream' );
        }
      }
      elsif ( $ch eq 'k3' ) {
        my $regex = "";
        $regex = $vs{'search_pattern'} if ( exists( $vs{'search_pattern'} ) );
        search( \%vs, $regex );
      }
      elsif ( $ch eq 'k2' ) {
        my $regex = "";
        $regex = $vs{'search_pattern'} if ( exists( $vs{'search_pattern'} ) );
        search( \%vs, $regex, -1 );
      }
      elsif ( $ch eq $ctrlg ) {
        my $newpos = $terminal->prompt_for_expression( \%vs, "Position?" );
        $vs{'cpos'} = localpos( \%vs, $newpos );
      }
      else{
        $refresh = 0;
      }
      $terminal->showVS( \%vs ) if ($refresh);
    }
  }
}

END{
  if ($CURRENT_INFO && -e $CURRENT_INFO){
    my $text = get_file_text($CURRENT_INFO);
    if ($text eq 'Locked'){
      print "Unlocking current gene...\n";
      open (OUT, ">$CURRENT_INFO") or die ("Could not unlock current gene\n");
      print OUT '';
      close OUT;
    }
  }
}

sub set_splice_frames {

  # Decides the frame for each exon/intron and intron/exon boundary
  # @_ = (hash_ref)
  # Returns nothing. Modifies the incoming hash by setting the correct sequence
  # and adding 'spliced_from' and 'spliced_to' tags (see 'modify')
  # For each exon/intron/exon junction, the procedure creates the sequence
  # derived from using each frame. The number of frames tried is determined by
  # the global variable $SPLICE_OVERHANG. Each of the resulting sequences is
  # blasted against the original query sequence. If the result of using a frame
  # yields the best blast score, this frame gets 7 points. If the result of
  # using a frame yields a 'gt' at the beginning of an intron or 'ag' at the
  # end of an intron, this frame gets 2 points. The frame with the most points
  # is chosen, and the sequence is modified accordingly.

  my $hash           = shift;
  my %points         = ();
  my @simils         = ();
  my @ids            = ();
  my @junction_frame = ();
  
  sort_exons($hash);
  my $n = $hash->{'features'}{'n'};
  
  ######### get gtag #########

  foreach my $exon ( keys %{ $hash->{'features'}{'exon'} } ) {
    if ( $hash->{'features'}{'exon'}{$exon}{'from'} > 0 ) {
      my $pre =
        $hash->{'features'}{'exon'}{$exon}{'from'} - ( $SPLICE_OVERHANG + 2 );
      my $post = $hash->{'features'}{'exon'}{$exon}{'to'};
      my $ag   = substr( $hash->{'seq'}, $pre, ( $SPLICE_OVERHANG + 2 ) );
      my $gt   = substr( $hash->{'seq'}, $post + 1, ( $SPLICE_OVERHANG + 2 ) );
      $hash->{'features'}{'exon'}{$exon}{'gt'} = $gt;
      $hash->{'features'}{'exon'}{$exon}{'ag'} = $ag;
    }
  }
  foreach my $exon ( keys %{ $hash->{'features'}{'exon'} } ) {
    if ( $exon > 0 ) {
      my $c = $SPLICE_OVERHANG;
      for ( my $i = -$c ; $i <= $c ; $i++ ) {
        my %totalseq = setseq( $hash, $exon, $i, $n );
        my $prot     = translate( $totalseq{'seq'} );
        my $fl       = gettempfile();
        open( OUT, ">$fl" );
        print OUT $prot;
        close OUT;
        my %temp = (
          'i' => $hash->{'aafile'},
          'j' => $fl,
          'p' => "blastp",
          'g' => 'T',
	        'G'  => '-1',
          'F' => 'F',
        );
        my %br = loadBlast( bl2seq( \%temp ) );
        $points{$i}{'ids'} = $br{0}{'hits'}{0}{'score'};
        push( @ids, $points{$i}{'ids'} );
        $points{$i}{'gt'} = $totalseq{'gt'};
        $points{$i}{'ag'} = $totalseq{'ag'};
        rename $fl, "trash";
      }
      my $mids = max(@ids);
      @ids = ();
#      debug($exon);
#      debug('-----');
      foreach my $junct ( keys %points ) {
        $points{$junct}{'total'} += 1 if ( $points{$junct}{'ids'} == $mids );
        $points{$junct}{'total'} += 3
          if ( lc( $points{$junct}{'gt'} ) eq "gt" );
        $points{$junct}{'total'} += 3
          if ( lc( $points{$junct}{'ag'} ) eq "ag" );
#        $points{$junct}{'total'} += 1 if ( $junct == 0 );
        $points{$junct}{'total'} = 0
          unless ( exists( $points{$junct}{'total'} ) );
#        debug($junct, $points{$junct}{'total'});
      }
      my @ordered =
        sort { $points{$b}{'total'} <=> $points{$a}{'total'} } keys %points;
      $junction_frame[$exon] = $ordered[0];
#      debug($ordered[0]);
#      debug('-----');
      %points = ();
    }
  }
  $hash->{'splice'} = \@junction_frame;
  modify( $hash, \@junction_frame, $n );
}

sub setseq {

  # Finds the sequence resulting from choosing a given splicing frame
  # @_ = ($hash_ref, exon_number, frame, total_number_of_exons)
  # Returns a hash containing the sequence, the two first bases of the
  # intron ('gt'), and the two last bases of the intron ('ag')

  my $hash   = shift;
  my $exon   = shift;
  my $i      = shift;
  my $n      = shift;
  my $temp = {};
  $temp = $hash->{'features'}{'exon'};
  my %result = ();
  for ( my $j = 0 ; $j < $exon - 1 ; $j++ ) {
    $result{'seq'} .= $temp->{$j}{'seq'} or die("$n\t$j");
  }
  my $e1 =
      $temp->{ $exon - 1 }{'seq'}
    . $temp->{ $exon - 1 }{'gt'};
  $e1 = $bftr->dright( $e1, $SPLICE_OVERHANG - $i );
  $result{'seq'} .= $bftr->dright( $e1, 2 );
  $result{'gt'} = $bftr->right( $e1, 2 );
  my $e2 =
      $temp->{$exon}{'ag'}
    . $temp->{$exon}{'seq'};
  $e2 = $bftr->dleft( $e2, $SPLICE_OVERHANG + $i );
  $result{'seq'} .= $bftr->dleft( $e2, 2 );
  $result{'ag'} = $bftr->left( $e2, 2 );

  for ( my $j = $exon + 1 ; $j <= $n ; $j++ ) {
    $result{'seq'} .= $temp->{$j}{'seq'};
  }
  return %result;
}

sub localize{
  # Set local blast
  my $hash    = shift;
  my $program = shift;
  my $suggest = shift;
  $suggest = 0 unless ($suggest);
  my $qfrom   = $hash->{'settings'}{'local'}{'contig'}{'tblastn'}{'query'}{'from'};
  my $qto     = $hash->{'settings'}{'local'}{'contig'}{'tblastn'}{'query'}{'to'};
  my $tfrom   = $hash->{'settings'}{'local'}{'contig'}{'tblastn'}{'template'}{'from'};
  my $tto     = $hash->{'settings'}{'local'}{'contig'}{'tblastn'}{'template'}{'to'};
  if ($qfrom && $qto && $tfrom && $tto &&
      $qto > $qfrom &&
      $tto > $tfrom)
  {
    foreach my $key (keys %{$hash->{'settings'}{$program}}){
      unless (exists($hash->{'settings'}{'local'}{$program}{$key})){
        $hash->{'settings'}{'local'}{$program}{$key} = $hash->{'settings'}{$program}{$key};
      }
    }
    $hash->{'settings'}{'local'}{$program}{'I'} = "$qfrom\,$qto";
    $hash->{'settings'}{'local'}{$program}{'J'} = "$tfrom\,$tto";
    if ($qto > $qfrom + 1){
      local_bl2seq($hash, $program, $suggest);
    }
    else{
      $terminal->message($hash, 'Query gap is not large enough');
    }
  }
  else{
    $terminal->message($hash, 'Invalid segment. No calculations performed.');
  }
}

sub local_bl2seq{
  # Perform local blast
  my $hash = shift;
  my $program = shift;
  my $suggest = shift;
  refreshVS($hash);
  my $init_e = $hash->{'settings'}{'local'}{$program}{'e'};
  my %hits = ();
  my $nhits = undef;
  %hits = loadBlast(bl2seq(\%{$hash->{'settings'}{'local'}{$program}}));
  $nhits = scalar keys %{$hits{0}{'hits'}};
#  $hash->{'settings'}{'local'}{$program}{'e'} *= 10 unless ($nhits);
  while ($hash->{'settings'}{'local'}{$program}{'e'} < 1e20 * $init_e && $suggest && !$nhits){
    %hits = loadBlast(bl2seq(\%{$hash->{'settings'}{'local'}{$program}}));
    $nhits = scalar keys %{$hits{0}{'hits'}};
    $hash->{'settings'}{'local'}{$program}{'e'} *= 10 unless ($nhits);
  }
  my $e = $hash->{'settings'}{'local'}{$program}{'e'};
  unless ($nhits){
    $nhits = 0;
    $hash->{'settings'}{'local'}{$program}{'e'} = $init_e;
  }
  $terminal->message($hash, "Local e-value: $e");
  $terminal->message($hash, "Added $nhits local hits");
  foreach my $hit (keys %{$hits{0}{'hits'}}){
    $bftr->hpush(\%{$hash->{'features'}{$program}{0}{'hits'}}, \%{$hits{0}{'hits'}{$hit}});
  }
}

sub gettempfile {

  # Returns a file name of the form 'tempn' (n=1,2,3,...) which
  # does not exist.

  my $result = "temp";
  my $n      = 0;
  while ( -e "$result$n" ) {
    $n++;
  }
  return "$result$n";
}




sub modify {

  #Modifies the sequence of the gene with the splicing info.

  my $hash   = shift;
  my $frames = shift;
  my $n      = shift;
  my $strand = $hash->{'strand'};
  for ( my $i = 1 ; $i <= $n ; $i++ ) {
    my $frame = $frames->[$i];
    my $j = $i - 1;
    $hash->{'features'}{'exon'}{$j}{'spliced_to'} =
      $hash->{'features'}{'exon'}{$j}{'to'} + $frame;
    if ( $frame < 0 ) {
      $hash->{'features'}{'exon'}{$j}{'seq'} =
        $bftr->dright( $hash->{'features'}{'exon'}{$j}{'seq'}, -1 * $frame );
    }
    elsif ( $frame > 0 ) {
      $hash->{'features'}{'exon'}{$j}{'seq'} .=
        $bftr->left( $hash->{'features'}{'exon'}{$j}{'gt'}, $frame );
    }
    $j = $i;
    $hash->{'features'}{'exon'}{$j}{'spliced_from'} =
      $hash->{'features'}{'exon'}{$j}{'from'} + $frame;
    if ( $frame > 0 ) {
      $hash->{'features'}{'exon'}{$j}{'seq'} =
        $bftr->dleft( $hash->{'features'}{'exon'}{$j}{'seq'}, $frame );
    }
    elsif ( $frame < 0 ) {
      $hash->{'features'}{'exon'}{$j}{'seq'} =
        $bftr->right( $hash->{'features'}{'exon'}{$j}{'ag'}, -1 * $frame )
        . $hash->{'features'}{'exon'}{$j}{'seq'};
    }
  }
}



#sub message{
#  my $hash = shift;
#  my $msg  = shift;
#  $bftr->hpush(\%{$hash->{'warnings'}}, $msg);
#}





######################################
##############Search##################
######################################

sub get_regex{
  #Prompts for regex and checks if valid.
  #Returns regex if valid or 0 otherwise.

  my $hash = shift;
  my $regex = $terminal->prompt_for_expression( $hash, "Pattern?" );
  #Check if regex is well formed
  my $check = 'anystring';
  eval '$check =~ /$regex/';
  if ($@){
    $terminal->message($hash, 'Incorrect regular expression');
    return 0;
  }
  $hash->{'search_pattern'} = $regex;
  return $regex;
}

sub search {

  # Looks for a pattern in the contig
  # @_ = (hash_ref, regex, direction)
  # hash_ref: \%vs
  # regex: string to search, in the form of a regular expression
  # direction: either 1 (forward) or -1 (backward)
  # Returns nothing. Points 'cpos' to the position of the hit

  my $hash  = shift;
  my $regex = shift;
  my $dir   = shift;
  $dir = 1 unless ($dir);
  my $cpos = $hash->{'cpos'};
  my $result = 0;
  if ( $dir != -1 ) {
    $result = getpos( \$hash->{'seq'}, $regex, $hash->{'cpos'} + 1 );
  }
  else {
    $result = getprevpos( \$hash->{'seq'}, $regex, $hash->{'cpos'} - 1 );
  }
  if ($result) {
    $hash->{'cpos'} = $result;
  }
  else {
    $terminal->message($hash, "Pattern not found");
  }
}

sub getpos {

  # Gets the position of the next pattern in a string
  # @_ = (string_ref, regex, start)
  # string_ref: reference to the template string
  # regex: expression to search (regular expression)
  # start: position of the template where the search starts
  # Returns the position of the first occurence of 'regex' in 'string'
  # after 'start'

  my $template = shift;    # By reference!
  my $target   = shift;
  my $ind      = shift;
  $ind = 0 unless ($ind);
  my $result = 0;
  pos($$template) = $ind;
  if ( $$template =~ /($target)/gi ) {
    $result = pos($$template) - length($1);
  }
  return $result;
}

sub getprevpos {

  # Gets the position of the previous pattern in a string
  # @_ = (string_ref, regex, start)
  # string_ref: reference to the template string
  # regex: expression to search (regular expression)
  # start: position of the template where the search starts
  # Returns the position of the first occurence of 'regex' in 'string'
  # before 'start'
  # This procedure follows the strategy 'Reverse the input! Reverse the regex!
  # Reverse the match!' by Jeff Pinyan (aka Japhy) in a simple way. Not fully
  # tested with complex expressions.

  my $template = shift;    # By reference!
  my $target   = shift;
  my $ind      = shift;
  $ind = 0 unless ($ind);
  my $rtemplate = substr( $$template, 0, $ind );
  $rtemplate = reverse($rtemplate);
  $target    = revregex($target);
  my $result = 0;
  pos($rtemplate) = 0;

  if ( $rtemplate =~ /$target/gi ) {
    $result = pos($rtemplate);
  }
  $result = $ind - $result if ($result);
  return $result;
}

sub backtranslate{
  # Takes a peptidic regex and returns the corresponding cDNA regex

  my $pept = shift;
  $pept = uc($pept);
  if ($pept =~ /(\[[^|\^]+?\])/){     # Or expressions: [xyz]
    my $choose = $1;                  # Changed to: (?:x|y|z)
    my $mod    = $choose;
    my $continue = 1;
    while ($continue){
      $continue = 0;
      $mod =~ s/([^\[\]])/\|$1/gi;  # OR between residues
      $mod =~ s/\[\|/\[/gi;         # no ORs at beginning or end
      $mod =~ s/\|\]/\]/gi;         #
      $mod =~ s/\[/\(\?\:/gi;       # Change parenthesis
      $mod =~ s/\]/\)/gi;           #
      $pept =~ s/\Q$choose\E/\Q$mod\E/i;
      $pept =~ s/\\//g;
      if ($pept =~ /(\[[^|]+?\])/){
        $choose = $1;
        $mod = $choose;
        $continue = 1;
      }
    }
  }
  my $result = '';
  my %bdict = (
    'A' => '(?:GC.)',
    'C' => '(?:TG[TC])',
    'D' => '(?:GA[TC])',
    'E' => '(?:GA[AG])',
    'F' => '(?:TT[TC])',
    'G' => '(?:GG.)',
    'H' => '(?:CA[TC])',
    'I' => '(?:AT[TCA])',
    'K' => '(?:AA[AG])',
    'L' => '(?:TT[AG]|CT.)',
    'M' => '(?:ATG)',
    'N' => '(?:AA[TC])',
    'P' => '(?:CC.)',
    'Q' => '(?:CA[AG])',
    'R' => '(?:AG[AG]|CG.)',
    'S' => '(?:AG[TC])',
    'T' => '(?:AC.)',
    'V' => '(?:GT.)',
    'W' => '(?:TGG)',
    'X' => '(?:.{3})',
    'Y' => '(?:TA[TC])',
    'Z' => '(?:TA[AG]|TGA)'
  );
  my @chars = split(//, $pept);
  foreach my $c (@chars){
    if (exists($bdict{$c})){
      $result .= $bdict{$c};
    }
    else{
      $result .= $c;
    }
  }
  return $result;
}

sub revregex {
  # Reverse the regex!
  my $regex  = shift;
  my $id     = shift;
  $id = '__1__' unless ($id);
  my @parts  = ();
  my @result = ();
  my $single_atom  = '[^\[\{\(]';
  my $escaped_atom = '\\+.';
  my $parenth_atom = '\([^()]*\)';
  my $other_atom   = '\[.*?\]|\{.*?\}';
  my $id_atom      = '__\d+__';
  if ($regex =~ /($parenth_atom)/){
    my $part = $1;
    my $continue = 1;
    while($continue){
      $continue = 0;
      $part = $bftr->dright($part, 1);
      $part = $bftr->dleft($part, 1);
      $regex =~ s/\(\Q$part\E\)/$id/;
      $id = getid($id);
      $part = revregex($part, $id);
      $part = "\($part\)";
      push (@parts, $part);
      if ($regex =~ /($parenth_atom)/){
        $part = $1;
        $continue = 1;
      }
    }
  }
  my $prefix = '';
  if ($regex =~ /^(.\:)/){
    $prefix = $1;
    $regex = $bftr->dleft($regex, length($1));
  }
  while ( $regex =~ /($escaped_atom|$other_atom|$id_atom|$single_atom)/g ) {
    my $atom = $1;
    my $mant = $bftr->left( $atom, 1 );
    if ( $mant eq '{' ||
         $mant eq '*' ||
         $mant eq '+' ) {
      my $patom = pop(@result);
      $atom = "$patom$atom";
    }
    if ( $bftr->left( $atom, 1) eq '('){
      $atom = $bftr->dleft($atom, 1);
      $atom = $bftr->dright($atom, 1);
      $atom = revregex($atom);
      $atom = "\($atom\)";
    }
    push( @result, $atom );
  }
  @result = reverse(@result);
  my $result = join( "", @result );
  if (@parts){
    my @rparts = reverse(@parts);
    $id = getid($id, 1);
    foreach my $part (@rparts){
      $result =~ s/$id/$part/;
      $id = getid($id, 1);
    }
  }
  $result = $prefix.$result if ($prefix);
  return $result;
}

sub getid{
  # Takes an ID from revregex (__N__) and returns the next (__N+1__) or
  # the previous (__N-1__) if $down=1
  my $pid  = shift;
  my $down = shift;
  my $result = '';
  $pid =~ s/_//g;
  if ($down){
    $result = $pid-1;
  }
  else{
    $result = $pid+1;
  }
  $result = "__$result\__";
  return $result;
}

######################################
###########Save result################
######################################

sub savegene {
  # Saves relevant information for further use
  my $hash = shift;
  my $xml  = '<?xml version="1.0" encoding="ISO-8859-1" ?>' . "\n";
  my $ts   = timestamp();
  my $info    = $bftr->path("$hash->{'fpath'}", "info.txt");
  my $hfile   = $bftr->path("$hash->{'fpath'}", "map.html");
  my $gfffile = $bftr->path("$hash->{'fpath'}", "$hash->{'file'}\.gff");
  $xml .= "<!-- Updated $ts -->\n";
  my $n    = $hash->{'features'}{'n'};
  my $cpos = 0;
  my $exon = '';
  my $gene = ">$hash->{'gname'}\n";
  my $save = {};
  my $wfile ='';
  # XML file
  for ( my $i = 0 ; $i <= $n ; $i++ ) {
    my $temp = $hash->{'features'}{'exon'}{$i};
    if ($i == 0){
      $temp->{'spliced_from'} = $temp->{'from'};
    }
    if ($i == $n){
      $temp->{'spliced_to'} = $temp->{'to'};
    }
    $save->{'gene'}{'exon'}{$i} = {
      'number' => $i + 1,
      'chromosome' => $hash->{'contig'}{'chr'},
      'template' => {
        'from' => realpos($hash, $temp->{'from'}),
        'to'   => realpos($hash, $temp->{'to'}),
        '_spliced_from' => realpos($hash, $temp->{'spliced_from'}),
        '_spliced_to'   => realpos($hash, $temp->{'spliced_to'})
      },
      'frame' => $hash->{'splice'}[$i],
    };
    if (exists($temp->{'warnings'})){
      my $c = $i + 1;
      $save->{'gene'}{'exon'}{$i}{'warnings'} = $temp->{'warnings'};
      my $w = $bftr->hjoin(', ', \%{$temp->{'warnings'}});
      $wfile .= "Exon $c\: $w\n";
    }
    $gene .= "$temp->{'seq'}\n";
  }
  $xml .= set_xml($save);
  #HTML map
  my $html = set_map($hash);
  #GFF file
  my $gff  = set_gff($hash);
  my $file = $hash->{'gfile'};
  open( OUT, ">$file" );
  print OUT $xml;
  close OUT;
  $file = $hash->{'sfile'};
  open( OUT, ">$file" );
  print OUT $gene;
  close OUT;
  $file = $hash->{'wfile'};
  open( OUT, ">$file" );
  print OUT $wfile;
  close OUT;
  open( OUT, ">$info" );
  print OUT "$ts\n";
  close OUT;
  open( OUT, ">$hfile" );
  print OUT $$html;
  close OUT;
  open( OUT, ">$gfffile" );
  print OUT $gff;
  close OUT;
}

################
##  GFF file  ##
################

sub set_gff{
  my $hash = shift;
  my $result = '';
  my $temp = \%{$hash->{'features'}{'exon'}};
  my $strand = '+';
  $strand = '-' if ($hash->{'strand'} eq '-1');
  foreach my $exon (sort{$temp->{$a}{'spliced_from'} <=> $temp->{$b}{'spliced_from'}} keys %{$temp}){
    my $f = realpos($hash, $temp->{$exon}{'spliced_from'});
    my $t = realpos($hash, $temp->{$exon}{'spliced_to'});
    my $from = min($f, $t);
    my $to   = max($f, $t);
    my $frame = 0;
    if (exists($hash->{'splice'}[$exon+1]) && $hash->{'splice'}[$exon+1]){
      $frame = $hash->{'splice'}[$exon+1] % 3;
    }
    $result .= "$hash->{'contig'}{'chr'}\t";
    $result .= "GeneTuner\t";
    $result .= "exon\t";
    $result .= "$from\t";
    $result .= "$to\t";
    $result .= "\.\t";
    $result .= "$strand\t";
    $result .= "$frame\t";
    $result .= "$hash->{'file'}\n";
  }
  return $result;
}


################
##  HTML map  ##
################

sub set_map{
  # Creates a map in HTML from hash info
  my $hash  = shift;
  my $template = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
<head>
  <title>__name__</title>
  <meta http-equiv="content-type" content="text/html; charset=iso-8859-1" />
<script language="JavaScript" type="text/javascript">
/*<![CDATA[*/
  function toggleIntron(intron){
    var n = document.getElementById(intron);
    if (n.className == "hidden"){
      n.className = "intron_collapse";
      var t = "toggle_" + intron;
      var toggle = document.getElementById(t);
      toggle.className = "hidden";
    }
    else{
      n.className = "hidden";
      var t = "toggle_" + intron;
      var toggle = document.getElementById(t);
      toggle.className = "intron_collapse";
    }
  }
  function showAll(){
    var dvs = document.getElementsByTagName(\'div\');
    for (var i = 0; i < dvs.length; i++){
      if (mleft(dvs[i].id, 5) == \'iExon\'){
        dvs[i].className = "intron_collapse";
        var t = "toggle_" + dvs[i].id;
        var toggle = document.getElementById(t);
        toggle.className = "hidden";
      }
    }
  }
  function hideAll(){
    var dvs = document.getElementsByTagName(\'div\');
    for (var i = 0; i < dvs.length; i++){
      if (mleft(dvs[i].id, 5) == \'iExon\'){
        dvs[i].className = "hidden";
        var t = "toggle_" + dvs[i].id;
        var toggle = document.getElementById(t);
        toggle.className = "intron_collapse";
      }
    }
  }
  function mleft(str, n){
	if (n <= 0)
	    return "";
	else if (n > String(str).length)
	    return str;
	else
	    return String(str).substring(0,n);
}

/*]]>*/
</script>
<style type="text/css">
  /*<![CDATA[*/
  a{
    cursor: pointer;
  }
  .text{
    font-family: Courier New;
    padding: 0;
    margin: 0;
    text-align: left;
    position: relative;
    float: left;
    text-indent: 0;
  }
  .exon{
    background-color: #CbCbCb;
    border: thin solid
  }
  .frame{
    background-color: #ffffa5;
    border: thin dotted;
  }
  .warning{
    font-family: Arial Narrow;
    background-color: #FAD2D2;
  	float: right;
    margin-left: 1em;
    border: thin outset;
  }
  .genome{
    font-weight: bold;
  }
  .rev{
    margin-bottom: 0.4em;
  }
  .transl{
  }
  .f3{
    margin-bottom: 1em;
  }
  .hidden{
    display: none;
  }
  .intron_collapse{
    cursor: pointer;
  }
  #header{
  	margin: 2em;
  	padding: 2em;
  }
  #map{
  	margin: 2em;
  	font-size: small;
  	padding:0;
  	width:90%;
  }
  /*]]>*/
</style>
</head>
<body>
  <div id="header">
    <span class="title">__title__</span>
  </div>
  <div id="map">
	  <div class="text">
	    __content__
	  </div>
	</div>
</body>
</html>
  ';
  my @estarts = ();
  my @eends   = ();
  my $debug = '';
  
  #Set titles
  my $chr = $hash->{'contig'}{'chr'};
  my $name = $hash->{'file'};
  my $ts = timestamp();
  my $strand = '+';
  $strand = '-' if ($hash->{'strand'} eq '-1');
  my $title = "<h1>Genomic map of $name<\/h1>".
              "<h2>Chromosome $chr \($strand strand\)<\/h2> Last updated: $ts\n<br \/><br \/>";
  foreach my $exon (sort{$a <=> $b} keys %{$hash->{'features'}{'exon'}}){
    my $temp = \%{$hash->{'features'}{'exon'}{$exon}};
    my $tfrom = $temp->{'spliced_from'};
    my $tto   = $temp->{'spliced_to'};
    push @estarts, $tfrom;
    push @eends, $tto;
    $temp->{'from'} = $tfrom;
    $temp->{'to'}   = $tto;
    my $e = $exon+1;
    $title .= '<a href="#Exon'.$e.'">Exon '.$e.'</a><br />'."\n";
  }
  $title .= '<a onclick="hideAll()">Hide long introns</a><br />'."\n";
  $title .= '<a onclick="showAll()">Show long introns</a><br />'."\n";
  my $start = min(@estarts) - $MAP_OVERHANG + int($MAP_WIDTH/2);
  my $end   = max(@eends)   + $MAP_OVERHANG - int($MAP_WIDTH/2);
  my $abs_start = realpos($hash, $start);
  my $abs_end   = realpos($hash, $end);
  my $nlen  = max(length($abs_start), length($abs_end));
  my $spacer = sp($nlen+1, ' ');
  my $nwindows = int(($end - $start)/$MAP_WIDTH) + 1;
  my $sep = sp($MAP_WIDTH, '-');
  $sep =~ s/(\-{9})\-/$1\+/g;
  #Set map
  my $map = '';
  my $exon_count = 0;
  my $orf_frame = 0;
  for (my $i = 0; $i <= $nwindows; $i++){
    $hash->{'cpos'} = $start + $i*$MAP_WIDTH;
    refreshS($hash, $MAP_WIDTH, 20, 1);
    my $gpos = realpos($hash, $start + $i*$MAP_WIDTH - int($MAP_WIDTH/2));
    $gpos = keep_length($nlen, $gpos);
    my $gline = $hash->{'screen'}{8};
    my $rev = greverse($gline);
    my @f = ();
    $f[2] = $hash->{'screen'}{12};
    $f[0] = $hash->{'screen'}{10};
    $f[1] = $hash->{'screen'}{11};
    my $warnings = '';
    if (exists($hash->{'highlights'})){  # Highlight exons
      my %hl = decodehl(\%{$hash->{'highlights'}});
      foreach my $k (keys %hl){    #Only highlights in the 8th line
        if ($hl{$k}{'y'} ne '8'){
          delete ($hl{$k});
        }
      }
      foreach my $h (sort{$hl{$b}{'x'} <=> $hl{$a}{'x'}} keys %hl){
        my $p = \%{$hash->{'features'}{'exon'}{$exon_count}};
        my $t = $hl{$h}{'text'};
        my $l = length($t);
        my $x = $hl{$h}{'x'};
        my $curr = $hash->{'cpos'} + $x - int($MAP_WIDTH/2);
        my $id = $exon_count + 1;
        $id = "Exon$id";
        my $next = $p->{'from'};
        if ($curr && $next && $curr == $next){
          $gline =~ s/(.{$x})(.{$l})/$1<span\_class\=\"exon\"_id\=\"$id\">$2<\/span>/;
          if (exists($p->{'warnings'})){
            $warnings = $bftr->hjoin('<br/>', \%{$p->{'warnings'}});
            $warnings = "<div_class=\"warning\">$warnings<\/div>";
          }
          $exon_count++ if (exists($hash->{'features'}{'exon'}{$exon_count+1}));
        }
        else{
          $gline =~ s/(.{$x})(.{$l})/$1<span\_class\=\"exon\">$2<\/span>/;
        }
        my $frame = ($curr) % 3;
        $frame = ($frame-$orf_frame) % 3;
        $f[$frame] =~ s/(.{$x})(.{$l})/$1<span\_class\=\"frame\">$2<\/span>/;
        $orf_frame = ($orf_frame+$l) % 3;
      }
      delete($hash->{'highlights'});
    }
    $map .= "<div_class=\"genome\">$warnings$spacer$gline<\/div>".
            "<div_class=\"sep\">$gpos\&nbsp;$sep<\/div>".
            "<div_class=\"genome_rev\">$spacer$rev<\/div>".
            "<div_class=\"transl_f1\">$spacer$f[2]<\/div>".
            "<div_class=\"transl_f2\">$spacer$f[0]<\/div>".
            "<div_class=\"transl_f3\">$spacer$f[1]<\/div>";
    $debug = '';
  }
  $map =~ s/\s/\&nbsp;/g;
  $map =~ s/<br\/>/<br \/>\n/g;
  $map =~ s/<\/div>/<\/div>\n/g;
  $map =~ s/\_/ /g;
  ###
  # Allow intron collapse
  my $left_limit  = '<span class="exon">';
  my $right_limit = '<span class="exon" id';
  my $llimit = '<div class="genome">';
  my $mline  = "$llimit.+?(?=$llimit)";

  my @introns = $bftr->isplit($right_limit, $map);
  $map = '';
  foreach my $intron (@introns){
    my @lines = $bftr->isplit($left_limit, $intron);
    if ($lines[0] =~ /<span class=\"exon\" id=\"(.+?)\"/){
      my $nintron = "i$1";
      my $nl = scalar @lines;
      my $hintron = '';
      if ($nl){
        $hintron = \$lines[$nl-1];
      }
      else{
        $hintron = \$lines[0];
      }
      my $len = 0;
      $len++ while ($$hintron =~ /$llimit/gsm);
      if ($len > $INTRON_COLLAPSE){   # Collapsible intron
        my $l = $len - 3;
        my $chunk = "((?:$mline))((?:$mline){$l})";
        my $coll_div = "<div class=\"intron_collapse\" id=\"$nintron\" ".
                       "onclick=\"toggleIntron('$nintron')\" title=\"Hide intron\">";
        $coll_div = "<div class=\"hidden\" id=\"toggle_$nintron\" ".
                    "onclick=\"toggleIntron('$nintron')\" title=\"Show intron\"><hr \/><\/div>".
                     $coll_div;
        $$hintron =~ s/$chunk/$1$coll_div$2<\/div>/gsm;
      }
    }
    $map .= join('', @lines);
  }
  ###
  
  #Place in template
  $template =~ s/__name__/$name/;
  $template =~ s/__title__/$title/;
  $template =~ s/__content__/$map/;
  return \$template;
}

sub greverse{
  # Get complementary sequence (not reversed)
  my $seq = shift;
  $seq =~ s/a/1/gi;
  $seq =~ s/c/2/gi;
  $seq =~ s/t/3/gi;
  $seq =~ s/g/4/gi;
  $seq =~ s/1/T/gi;
  $seq =~ s/2/G/gi;
  $seq =~ s/3/A/gi;
  $seq =~ s/4/C/gi;
  return $seq;
}

sub keep_length{
  # Right-justify with spaces a number $n2 so that its length
  # is $n1.
  my $n1 = shift;
  my $n2 = shift;
  my $result = $n2;
  if (length($n2) < $n1){
    my $l = $n1 - length($n2);
    my $s = sp($l, ' ');
    $result = $s.$n2;
  }
  return $result;
}

####################
##  end HTML map  ##
####################

sub timestamp {

# Adapted from Kirk Brown: "Using localtime to find the time in your Perl scripts"
# http://perl.about.com/od/perltutorials/a/perllocaltime.htm

  my @months   = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
  my (
    $second,     $minute,    $hour,
    $dayOfMonth, $month,     $yearOffset,
    $dayOfWeek,  $dayOfYear, $daylightSavings
    )
    = localtime();
  my $year    = 1900 + $yearOffset;
  my $theTime =
"$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
  return $theTime;
}

sub sort_exons {
  # Sorts exons
  my $hash = shift;
  my $temp = {};
  my $n = 0;
  foreach my $sorted ( sort {
  $hash->{'features'}{'exon'}{$a}{'from'} <=> $hash->{'features'}{'exon'}{$b}{'from'}
  } keys %{ $hash->{'features'}{'exon'} }) {
    $n = $bftr->hpush($temp, $hash->{'features'}{'exon'}{$sorted});
  }
  $hash->{'features'}{'exon'} = $temp;
  $hash->{'features'}{'n'} = $n - 1;
  
}

sub setundo {
  # Sets a copy of 'exon' subhash to allow its recovery
  my $hash = shift;
  %{ $hash->{'undo'} } =
    %{ Storable::dclone( \%{ $hash->{'features'}{'exon'} } ) };
}

sub nextfeat {
  #
  my $hash   = shift;
  my $cpos   = shift;
  my $result = -1;
  foreach my $n (%$hash) {
    my $to = $hash->{$n}{'to'} if exists( $hash->{$n} );
    if ( $to && $to > $cpos ) {
      if ( $result == -1 ) {
        $result = $n;
      }
      else {
        $result = $n if ( $to < $hash->{$result}{'to'} );
      }
    }
  }
  return $result;
}

sub prevfeat {
  my $hash   = shift;
  my $cpos   = shift;
  my $result = -1;
  foreach my $n (%$hash) {
    my $from = $hash->{$n}{'from'} if exists( $hash->{$n} );
    if ( $from && $from < $cpos ) {
      if ( $result == -1 ) {
        $result = $n;
      }
      else {
        $result = $n if ( $from > $hash->{$result}{'from'} );
      }
    }
  }
  return $result;
}

sub getprevM {
  my $hash   = shift;
  my $i      = $hash->{'cpos'};
  my $result = rindex( lc( $hash->{'seq'} ), 'atg', $i - 1 );
  return $result;
}

sub getcexon {
  my $hash = shift;
  my $p    = $hash->{'cpos'};
  $hash->{'current_exon'} = -1;
  foreach my $n ( keys %{ $hash->{'features'}{'exon'} } ) {
    $hash->{'current_exon'} = $n
      if overlap( $hash->{'features'}{'exon'}{$n}, $p );
  }
}

sub delexon {

  # Delete exon

  my $hash = shift;
  getcexon($hash);
  my $c = $hash->{'current_exon'};
  delete $hash->{'features'}{'exon'}{$c}
    if exists( $hash->{'features'}{'exon'}{$c} );
}

sub insertexon {

 # Starts an exon at the current position or edits an existing one
 # @_ = (hash_ref)
 # Returns nothing
 # If 'cpos' is inside a previously existing exon, this exon is edited, with
 # the same starting position. Else, a new exon is created. Exons are not sorted
 # in any way. The sorting must be carried out before saving.

  my $hash = shift;
  return if ( $hash->{'step'} == 3 );
  $hash->{'step'} = 3;
  my $n = 0;
  foreach my $k ( keys %{ $hash->{'features'}{'exon'} } ) {
    if ( overlap( \%{ $hash->{'features'}{'exon'}{$k} }, $hash->{'cpos'} ) ) {
      my $cpos = $hash->{'features'}{'exon'}{$k}{'from'};
      $hash->{'cpos'} = $hash->{'features'}{'exon'}{$k}{'from'} + 3
        if ( $hash->{'cpos'} <= $cpos + 3 );
      my $mod = ( $hash->{'cpos'} - $cpos ) % 3;
      my $sum = 3 - $mod if ($mod);
      $hash->{'cpos'} += $sum if ($sum);
      $hash->{'features'}{'exon'}{$k}{'to'} = $hash->{'cpos'} - 1;
      $hash->{'current_exon'} = $k;
      return;
    }
  }
  while ( exists( $hash->{'features'}{'exon'}{$n} ) ) {
    $n++;
  }
  $hash->{'features'}{'exon'}{$n}{'from'} = $hash->{'cpos'};
  $hash->{'cpos'} += 3;
  $hash->{'features'}{'exon'}{$n}{'to'} = $hash->{'cpos'} - 1;
  $hash->{'current_exon'} = $n;
}

sub getnext {

  # Gets the position of the next feature

  my $h       = shift;
  my $cpos    = shift;
  my $result  = $cpos;
  my @numbers = ();
  foreach my $n ( keys %$h ) {
    my $diff = $h->{$n}{'from'} - $cpos;
    my $dift = $h->{$n}{'to'} - $cpos;
    push @numbers, $diff if ( $diff > 0 );
    push @numbers, $dift if ( $dift > 0 );
  }
  my @o = sort { $a <=> $b } (@numbers);
  $result += $o[0] if ( $o[0] );
  return $result;
}

sub getprev {

  # Gets the position of the previous feature

  my $h       = shift;
  my $cpos    = shift;
  my $result  = $cpos;
  my @numbers = ();
  foreach my $n ( keys %$h ) {
    my $diff = $h->{$n}{'from'} - $cpos;
    my $dift = $h->{$n}{'to'} - $cpos;
    push @numbers, $diff if ( $diff < 0 );
    push @numbers, $dift if ( $dift < 0 );
  }
  my @o = sort { $b <=> $a } (@numbers);
  $result += $o[0] if ( $o[0] );
  return $result;
}

sub initVS {

#VS->{'settings'}->{'tblastn'}->...          ----> tblasn settings (-i, -j, ...)
#                ->{'blastn'} ->...          ----> blasn settings (-i, -j, ...)
#  ->{'features'}->{'exon'}  ->{n}->{'from'}
#                                  ->{'to'}
#               ->{'tblastn'}->loadBlast...  ----> see loadBlast
#               ->{'blastn'} ->loadBlast...
#  ->{'cpos'}                                ----> current position
#  ->{'contig'}                              ----> contig start and end
#  ->{'seq'}                                 ----> contig sequence

  my $hash  = shift;
  my $file  = $hash->{'file'};
  my $ppath = $settings{'ppath'};
  my $pext  = $settings{'aa_ext'};
  my $npath = $settings{'npath'} if (exists($settings{'npath'}));
  my $ntext = $settings{'nt_ext'} if (exists($settings{'nt_ext'}));
  my $tpath = $settings{'dbpath'};
  print "Tying database\n";
  print "This may take a few minutes the first time...\n";
  my $tobj = tie %db, 'Bio::DB::Fasta', $tpath;
  if (exists($hash->{'contig'})){
    foreach my $n (keys %{$hash->{'features'}{'exon'}}){
      $hash->{'cpos'} = $hash->{'features'}{'exon'}{$n}{'from'};
      $hash->{'features'}{'exon'}{$n}{'from'} = realpos($hash);
      $hash->{'cpos'} = $hash->{'features'}{'exon'}{$n}{'to'};
      $hash->{'features'}{'exon'}{$n}{'to'}   = realpos($hash);
    }
  }
  else{
    get_gene($hash);
  }
  getfromto( $hash, $CONTIG_OVERHANG );
  $hash->{'strand'} =
    ( $hash->{'contig'}{'from'} < $hash->{'contig'}{'to'} ) ? 1 : -1;
  sort_exons($hash);
  my $seq  = '';
  my $tlen = tied(%db)->length( $hash->{'contig'}{'chr'} );
  $hash->{'contig'}{'from'} = 1     if ( $hash->{'contig'}{'from'} < 1 );
  $hash->{'contig'}{'from'} = $tlen if ( $hash->{'contig'}{'from'} > $tlen );
  $hash->{'contig'}{'to'}   = 1     if ( $hash->{'contig'}{'to'} < 1 );
  $hash->{'contig'}{'to'}   = $tlen if ( $hash->{'contig'}{'from'} > $tlen );

#  while ( !$seq ) {
    $seq =
      $db{
"$hash->{'contig'}{'chr'}:$hash->{'contig'}{'from'},$hash->{'contig'}{'to'}"
      }
      or die(
"$file\t$hash->{'contig'}{'chr'}:$hash->{'contig'}{'from'},$hash->{'contig'}{'to'}"
      );
#  }

  my $fname = $bftr->strip_file_name($file);
  my $from  = min( $hash->{'contig'}{'from'}, $hash->{'contig'}{'to'} );
  my $to    = max( $hash->{'contig'}{'from'}, $hash->{'contig'}{'to'} );
  $hash->{'contig'}{'from'} = $from;
  $hash->{'contig'}{'to'}   = $to;
  my $p = $bftr->path($settings{'basepath'}, $fname, "$fname\_contig.seq");
  $hash->{'gname'} = $file;
  $hash->{'db'}    = $tpath;
  open( OUT, ">$p" );
  print OUT ">$file\n$seq";
  close OUT;
  my $aafile = $bftr->path($ppath, "$fname.$pext");
  $aafile =~ s/\.\./\./g;
  my $ntfile = $bftr->path($npath, "$fname.$ntext");
  $ntfile =~ s/\.\./\./g;
  my $qu     = Bio::SeqIO->new( -file => $aafile, '-format' => 'Fasta' );
  my $qobj   = $qu->next_seq();
  my $l      = $qobj->length();
  $hash->{'qlength'} = $l;
  $hash->{'aafile'}  = $aafile;
  my %tdefault = (
    'i' => $aafile,
    'j' => $p,
    'p' => "tblastn",
    'F' => 'F',
    'g' => 'T',
##    'G' => -1,
#    'E' => -1,
##    'W' => 3,
    'e' => 10
  );
  my %ndefault = (
    'i' => $ntfile,
    'j' => $p,
    'p' => "blastn",
    'g' => 'T',
#    'G' => -1,
#    'E' => -1,
#    'W' => 11,
    'F' => 'F',
    'e' => 10
  );
  %{ $hash->{'settings'}->{'tblastn'} } = %tdefault;
  %{ $hash->{'settings'}->{'blastn'} }  = %ndefault;
  $hash->{'settings'}{'fcol'} = 0;
  refreshVS($hash);
  refreshgene($hash);
  $hash->{'seq'} = $seq;
  my @g = keys( %{ $hash->{'features'}{'exon'} } );
  my $f = min(@g);
  if ($hash->{'strand'} == -1){
    $f = max(@g);
  }
  $hash->{'cpos'} = $hash->{'features'}{'exon'}{$f}{'from'};
  $hash->{'step'} = 1;
}

sub refreshgene {
  my $h = shift;
#  debug(hash_to_string(\%{$h->{'features'}{'exon'}}));
  foreach my $e (keys %{$h->{'features'}{'exon'}}){
    $h->{'features'}{'exon'}{$e}{'from'} = localpos($h, $h->{'features'}{'exon'}{$e}{'from'});
    $h->{'features'}{'exon'}{$e}{'to'}   = localpos($h, $h->{'features'}{'exon'}{$e}{'to'});
  }
}

sub refreshVS {
  my $hash = shift;
  ########### Refresh blasts ############
  my $tblastn = bl2seq( \%{ $hash->{'settings'}->{'tblastn'} } );
  my $blastn  = bl2seq( \%{ $hash->{'settings'}->{'blastn'} } );
  %{ $hash->{'features'}->{'tblastn'} } = loadBlast($tblastn);
  %{ $hash->{'features'}->{'blastn'} }  = loadBlast($blastn);
}

sub clearVS {
  my $hash = shift;
  $hash->{'screen'} = ();
}

sub refreshS {
  # Fills the 'screen' subhash
  my $hash   = shift;
  my $ncols  = shift;
  my $nrows  = shift;
  my $nopad  = shift;
  my $rrow   = int( $nrows / 2 );
  my $pad    = int( $ncols / 12 );
  $pad = 0      if ($nopad);
  my $win    = int( 5 * $ncols / 6 );
  $win = $ncols if ($nopad);
  my $center = int( $ncols / 2 );
  my $line   = pad($ncols);
  clearVS($hash);
  ########### Limit position ###########
  my $from = min( $hash->{'contig'}{'from'}, $hash->{'contig'}{'to'} );
  my $to   = max( $hash->{'contig'}{'from'}, $hash->{'contig'}{'to'} );
  if ( $hash->{'cpos'} < 1 ) {
    $hash->{'cpos'}     = 1;
    $terminal->message($hash, "Reached beginning");
  }
  if ( $hash->{'cpos'} > $to - $from ) {
    $hash->{'cpos'}     = $to - $from;
    $terminal->message($hash, "Reached end");
  }
  ###########   Headers    ###########
  my $query_name   = $hash->{'file'};
  my $query_length = $hash->{'qlength'};
  my $contig_name  = $hash->{'contig'}{'chr'};
  my $strand       = '+';
  $strand = '-' if ($hash->{'strand'} eq '-1');
  
  ########### Current exon ###########
  if ( $hash->{'step'} == 3 ) {    ########### Exon mode ###########
    my $cexon = $hash->{'current_exon'};
    my $cpos  = $hash->{'features'}{'exon'}{$cexon}{'from'};
    if ( $hash->{'cpos'} < ( $hash->{'features'}{'exon'}{$cexon}{'from'} + 3 ) )
    {
      $hash->{'cpos'} = $hash->{'features'}{'exon'}{$cexon}{'from'} + 3;
    }
    my $mod = ( $hash->{'cpos'} - $cpos ) % 3;
    my $sum = 3 - $mod if ($mod);
    $hash->{'cpos'} += $sum if ($sum);
    $hash->{'features'}{'exon'}{$cexon}{'to'} = $hash->{'cpos'} - 1;
  }
  else {                           ########### Normal mode ###########
    getcexon($hash);
    my $cx = $hash->{'current_exon'};
    if ( $cx > -1 && exists($hash->{'features'}{'exon'}{$cx}{'warnings'}) ) {
      foreach my $w (sort{$a <=> $b} keys %{$hash->{'features'}{'exon'}{$cx}{'warnings'}}){
        $terminal->message($hash, $hash->{'features'}{'exon'}{$cx}{'warnings'}{$w});
      }
    }
  }
  ########### Load exons ###########
  foreach my $nex ( keys %{ $hash->{'features'}{'exon'} } ) {
    my $f  = $hash->{'features'}{'exon'}{$nex}{'from'};
    my $t  = $hash->{'features'}{'exon'}{$nex}{'to'};
    $hash->{'features'}{'exon'}{$nex}{'qtext'} =
      subseq( \$hash->{'seq'}, $f, $t - $f + 1 )
      or die length(\$hash->{'seq'}).": $f\t$t\n";
    $hash->{'features'}{'exon'}{$nex}{'seq'} =
      $hash->{'features'}{'exon'}{$nex}{'qtext'};
  }
  my $f = $hash->{'cpos'} - $center + $pad;
  ########### Load sequence ############
  my $ttr = subseq( \$hash->{'seq'}, $f - 2, $win + 4 );

  #my $seq = substr( $hash->{'seq'}, $f, $win );

  my $seq = subseq( \$ttr, 2, $win );
  my %screen = ();
  $screen{'from'} = $f;
  $screen{'to'}   = $f + $win;
  ########### Translations ############
  $ttr = subseq( \$hash->{'seq'}, $f - 2, $win + 4 );
  my $frame = ( $f + 1 ) % 3;
  my @frames = get_frames( $ttr, $frame );
  ########### Load blastn features ############
  my $showfeat =
    getfeats( \%{ $hash->{'features'}{'blastn'}{0}{'hits'} }, \%screen );
  my @bnfeat = loadfeats( $showfeat, \%screen );
  my $bnt = sticks( $seq, $bnfeat[0] );
  ########### Load tblastn features ############
  $showfeat = {};
  $showfeat =
    getfeats( \%{ $hash->{'features'}{'tblastn'}{0}{'hits'} }, \%screen );
  my @tbnfeat = loadfeats( $showfeat, \%screen, 1 );

  #$tbnfeat = insert($line, $pad, $tbnfeat);
  my @tbnt = ();
  my $tsticks = $line;
  foreach my $tfeat (@tbnfeat){
    $tbnt[0] = sticks( $frames[0], $tfeat );
    $tbnt[1] = sticks( $frames[1], $tfeat );
    $tbnt[2] = sticks( $frames[2], $tfeat );
    foreach my $l (@tbnt) {
      $tsticks = insert( $tsticks, $pad, $l, ' ' ) if ( $l =~ /\|/ );
    }
  }
  ########### Position in the chromosome ############
  my $p = realpos($hash);
  ########### Load highlights ############
  # Exons
  $showfeat = {};
  $showfeat = getfeats( \%{ $hash->{'features'}{'exon'} }, \%screen );
  my @exs = loadfeats( $showfeat, \%screen );
  $exs[0] = insert( $line, $pad, $exs[0] );
  my $n = $bftr->getnkeys( \%{ $hash->{'highlights'} } );
  $hash->{'highlights'}{$n}{'text'} = $exs[0];
  $hash->{'highlights'}{$n}{'line'} = $rrow - 2;

  # Current residue
  $n = $bftr->getnkeys( \%{ $hash->{'highlights'} } );
  my $hf           = $hash->{'cpos'} % 3;
  my $currentcodon = subseq( \$hash->{'seq'}, $hash->{'cpos'}, 3 );
  my $currentres   = translate($currentcodon);
  $hash->{'highlights'}{$n}{'text'} = insert( $line, $center, $currentres );
  $hash->{'highlights'}{$n}{'line'} = $rrow + $hf;
  $hash->{'highlights'}{$n}{'effect'} = 2;    #reverse, bold
  ########### Refresh screen ############
  $hash->{'screen'}->{0} = "$query_name\tlength\: $query_length";
  $hash->{'screen'}->{1} = "Chromosome: $contig_name\tstrand\: $strand";
  $hash->{'screen'}->{ $rrow - 6 } = insert( $line, $center, "\|$p" );
  $hash->{'screen'}->{ $rrow - 5 } = insert( $line, $pad,    $bnfeat[1] );
  $hash->{'screen'}->{ $rrow - 4 } = insert( $line, $pad,    $bnfeat[0] );
  $hash->{'screen'}->{ $rrow - 3 } = insert( $line, $pad,    $bnt );
  $hash->{'screen'}->{ $rrow - 2 } = insert( $line, $pad,    $seq );
  $hash->{'screen'}->{ $rrow + 0 } = insert( $line, $pad,    $frames[0] );
  $hash->{'screen'}->{ $rrow + 1 } = insert( $line, $pad,    $frames[1] );
  $hash->{'screen'}->{ $rrow + 2 } = insert( $line, $pad,    $frames[2] );
  $hash->{'screen'}->{ $rrow + 3 } = $tsticks;
  $hash->{'screen'}->{ $rrow + 4 } = insert( $line, $pad,    $tbnfeat[0] );
  $hash->{'screen'}->{ $rrow + 5 } = insert( $line, $pad,    $tbnfeat[1] );
  $hash->{'screen'}->{ $rrow + 6 } = insert( $line, $pad,    $tbnfeat[2] ) if ($tbnfeat[2]);
  $hash->{'screen'}->{ $rrow + 7 } = insert( $line, $pad,    $tbnfeat[3] ) if ($tbnfeat[3]);

  if ( exists( $hash->{'warnings'} ) ) {
    $hash->{'screen'}->{ $nrows - 2 } = $bftr->hjoin('. ', \%{$hash->{'warnings'}});
    delete $hash->{'warnings'};
  }
  else {
    $hash->{'screen'}->{ $nrows - 2 } = '';
  }
}

sub subseq {

  # Gets a part of a sequence, padding with Ns if necessary
  # @_ = (seq_ref, pos, length)
  # Returns a substring of seq starting at pos with length
  # nucleotides.
  # pos and length can have any numerical value. If they are
  # outside seq, the rest is padded with Ns

  my $seq    = shift;
  my $pos    = shift;
  my $len    = shift;
  my $totl   = length($$seq);
  my $result = '';
  if ( $pos < 0 && $pos + $len > $totl ) {
    $result = pad( -$pos, 'N' );
    $result .= $$seq;
    $result .= pad( $pos + $len - $totl, 'N' );
  }
  elsif ( $pos < 0 ) {
    $result = pad( -$pos, 'N' );
    $result .= substr( $$seq, 0, $pos + $len );
  }
  elsif ( $pos + $len > $totl ) {
    my $l = $totl - $pos;
    if ($pos<0 || $l>length($$seq) || $l < 0){
      die("$$seq\t$pos\t$totl\n");
    }
    $result = substr( $$seq, $pos, $totl - $pos );
    $result .= pad( $pos + $len - $totl, 'N' );
  }
  else {
    $result = substr( $$seq, $pos, $len );
  }
  return $result;
}

sub realpos {
  my $hash = $_[0];
  my $cpos = '';
  if ( $_[1] ) {
    $cpos = $_[1];
  }
  else {
    $cpos = $hash->{'cpos'};
  }
  my $result = 0;
  if ( $hash->{'strand'} == 1 ) {
    $result = $hash->{'contig'}{'from'} + $cpos;
  }
  else {
    $result = $hash->{'contig'}{'to'} - $cpos;
  }
  return $result;
}

sub localpos {
  my $hash = shift;
  my $cpos = shift;
  $cpos = $hash->{'cpos'} unless ($cpos);
  my $result = 0;
  if ( $hash->{'strand'} == 1 ) {
    $result = $cpos - $hash->{'contig'}{'from'};
  }
  else {
    $result = $hash->{'contig'}{'to'} - $cpos;
  }
  return $result;
}

sub get_frames {
  my @result = '';
  my $seq    = shift;
  my $frame  = shift;
  my $l      = length($seq);
  my $win    = $l - 4;
  for ( my $i = 0 ; $i < 3 ; $i++ ) {
    my $n    = ( $frame + $i ) % 3;
    my $sq   = substr( $seq, $i, $l - $i );
    my $line = translate($sq);
    $line       = pad( $i + 1 ) . pformat($line);
    $line       = substr( $line, 3, $win );
    $result[$n] = $line;
  }
  return @result;
}

sub getfeats {

  # Gets the gene features which are inside the screen
  # @_ = (exons_hash_ref, screen_hash_ref)
  # Returns a hash containing each gene feature which
  # overlaps with the screen

  my $h1     = shift;
  my $h2     = shift;
  my $result = {};
  foreach my $n ( keys %$h1 ) {
    my $t = $bftr->getnkeys( $result );
    $result->{$t} = $h1->{$n} if ( overlap( \%{ $h1->{$n} }, $h2 ) );
  }
  return $result;
}

sub sticks {

# Creates lines with sticks between aligned sequences
# @_ = (line1, line2)
# Returns a string with sticks (|) at the positions where <line1> equals <line2>

  my $t1     = shift;
  my $t2     = shift;
  my @result = ();
  my $result = '';
  my @t1     = split( //, $t1 );
  my @t2     = split( //, $t2 );
  my $i1     = @t1;
  my $i2     = @t2;
  my $t      = min( $i1, $i2 );
  $result = pad($t);

  for ( my $i = 0 ; $i < $t ; $i++ ) {
    if ( $t1[$i] =~ /\w/ ) {
      $result = insert( $result, $i, "|" )
        if ( uc( $t1[$i] ) eq uc( $t2[$i] ) );
    }
  }
  return $result;
}

sub min {
  my @f      = @_;
  my $result = $f[0];
  foreach my $f2 (@f) {
    my $f = $f2;
    $f2 =~ s/\D\-//g;
    if ( $f eq $f2 ) {
      $result = $f if ( $f < $result );
    }
  }
  return $result;
}

sub max {
  my @f      = @_;
  my $result = $f[0];
  foreach my $f2 (@f) {
    my $f = $f2;
    $f2 =~ s/\D\-//g;
    if ( $f eq $f2 ) {
      $result = $f if ( $f > $result );
    }
  }
  return $result;
}

#########################################
sub debug {
  my @deb = @_;
  open( OUT, ">>debug.txt" );
  foreach my $deb (@deb) {
    print OUT "$deb\t";
  }
  print OUT "\n";
  close OUT;
}

sub phash {
  my $h      = $_[0];
  my $result = '';
  foreach my $key ( keys %$h ) {
    $result .= "$key\t$h->{$key}\n";
  }
  return $result;
}
#########################################

sub insert {

 # Inserts a text inside another text at a specific position
 # @_ = (canvas, pos, text, transp)
 # Returns the <canvas> after overprinting <text> characters except <transp>
 # starting from <pos> and trimming any overflow.
 # Useful for displaying data in a line. Usually, the line is started with
 # a canvas made of blank spaces with the 'sp' function and filled with the data
 # to display. The resulting string can be used as a canvas to diplay further
 # information.

  my $canvas = shift;
  my $cpos   = shift;
  my $text   = shift;
  my $transp = shift;
  $transp = '' unless($transp);
  my $result = $canvas;
  my $l      = length($text);
  if ( $cpos < -$l || $cpos > length($canvas) ) {
    return $canvas;
  }
  if ( $cpos + $l > length($canvas) ) {
    $l = length($canvas) - $cpos;
    $text = substr( $text, 0, $l );
  }
  if ( $cpos < 0 ) {
    $l += $cpos;
    $text = substr( $text, -$cpos, $l );
    $cpos = 0;
  }
  if ($transp){
    for (my $p = 0; $p < $l; $p++){
      my $c = substr $text, $p, 1;
      unless ($c eq $transp){
        substr($result, $cpos + $p, 1, $c);
      }
    }
  }
  else{
    substr $result, $cpos, $l, $text;
  }
  return $result;
}

sub loadfeats {

  # Inserts a feature into a line
  # @_ = (feat_hashref, screen_hashref, isprot)
  # Returns a string ready to display
  # <feat> contains a number of features, with 'from', 'to', and 'qtext' entries
  # <screen> is a subhash of VS
  # <isprot> is used as a boolean. If true, 'qtext' is considered a protein,
  # and spaces are inserted to align with DNA lines.

  my $features = shift;
  my $screen   = shift;
  my $isprot   = shift;
  my $l        = $screen->{'to'} - $screen->{'from'};
  my $line     = pad($l);
  my @result = ();
  $result[0]   = $line;
  $result[1]   = $line;
  my @k = keys %{$features};
  my $keys = \@k;
  if (@k && exists($features->{$k[0]}->{'expect'})){
    my @temp = sort{$features->{$b}{'score'} <=> $features->{$a}{'score'}} @$keys;
    $keys = \@temp;
  }
  my @prev = ();
  foreach my $n ( @$keys ) {
    my $offset = $features->{$n}->{'from'} - $screen->{'from'};
    my $pl     = $features->{$n}->{'to'} - $features->{$n}->{'from'};
    my $text   = $features->{$n}->{'qtext'};
    my $end    = $offset + length($text);
    my $add    = 0;
    if (is_there_overlap(@prev, \%{$features->{$n}})){
      $add = 2;
    }
    if ($isprot){
      $end  = $offset + 3 * (length($text) - 1);
      $text = pformat($text);
    }
    $result[0 + $add]   = $line unless ($result[0 + $add]);
    $result[1 + $add]   = $line unless ($result[1 + $add]);
    $result[0 + $add] = insert( $result[0 + $add], $offset, $text );
    # Numbers in query
    if (exists($features->{$n}{'query'}{'from'})){
      my $from = $features->{$n}{'query'}{'from'};
      my $to   = $features->{$n}{'query'}{'to'};
      $result[1 + $add] = insert($result[1 + $add], $offset, $from);
      $result[1 + $add] = insert($result[1 + $add], $end   , $to);
      push @prev, \%{$features->{$n}} unless ($add);
    }
  }
  return @result;
}

sub is_there_overlap{

  # Checks overlaps in a number of hashes and subhashes
  # @_ = (list_of_hash_refs)
  # Returns 1 if there is an overlap of the last hash_ref with any of the previous ones,
  # 0 otherwise
  # Each hash_ref must point to a hash containing keys named 'from' and 'to'
  
  my @hashes = @_;
  my @keys = ();
  my $result = '';
  my $h2 = $hashes[@hashes-1];
  if (exists($h2->{'from'}) && exists($h2->{'to'})){
    for (my $i = 0; $i < @hashes - 1; $i++){
      my $h1 = $hashes[$i];
      if (exists($h1->{'from'}) && exists($h1->{'to'})){
        if (overlap($h1, $h2)){
          $result = 1;
        }
      }
    }
  }
  return $result;
}

sub delinsert {

  # Deletes bases in <query> if they align to gaps in <template>

  my $query    = shift;
  my $template = shift;
  my @qu       = split( //, $query );
  my @te       = split( //, $template );
  my $i1       = @qu;
  my $i2       = @te;
  my $t        = min( $i1, $i2 );
  my $result   = '';
  for ( my $i = 0 ; $i < $t ; $i++ ) {
    $result .= $qu[$i] if ( $te[$i] ne '-' );
  }
  return $result;
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

sub pad {

  # Produces a string consisting of a given number of spaces
  # @_ = (length_of_string, padding_char)
  # Returns a scalar with <length_of_string> padding chars
  # By default, <padding_char> is a space.

  my $n = shift;
  my $c = shift;
  $c = ' ' unless ($c);
  return '' unless ($n);
  my $result = '';
  for ( my $i = 0 ; $i < $n ; $i++ ) {
    $result .= $c;
  }
  return $result;
}

sub decodehl {
  my $hl     = shift;
  my %result = ();
  my $eff    = 0;
  foreach my $n ( keys %$hl ) {
    if ( exists( $hl->{$n}{'effect'} ) ) {
      $eff = $hl->{$n}{'effect'};
    }
    else {
      $eff = 0;
    }
    my $y = $hl->{$n}{'line'};
    my $t = $hl->{$n}{'text'};
    while ( $t =~ /([^\s]+)/g ) {
      my $text   = $1;
      my $offset = pos($t) - length($text);
      my $n      = $bftr->getnkeys( \%result );
      $result{$n}{'x'}      = $offset;
      $result{$n}{'y'}      = $y;
      $result{$n}{'text'}   = $text;
      $result{$n}{'effect'} = $eff;
    }
  }
  return %result;
}

sub bl2seq {

  # Executes the bl2seq program with options
  # @_ = (hash_ref)
  # Returns a scalar with the text piped from the execution
  # The routine reads the options from the hash directly. Hash_ref
  # points to the hash in the 'settings' hash that contains blast
  # options

  my $hash    = shift;
  my $command = "bl2seq ";
  foreach my $argument ( keys %$hash ) {
    $command .= " \-$argument $hash->{$argument}";
  }
  my $result = `$command`;
  return $result;
}

sub get_gene {

  # Gets all of the previously stored information about the gene
  # @_ = ($hash_ref)
  # Returns no value
  # Modifies the passed hash to include the
  # {'features'}{'exon'} subhash, plus
  # 'fpath', 'gfile', 'sfile', and 'wfile'

  my $hash     = shift;
  my $file     = $bftr->path("$settings{'basepath'}", "$hash->{'file'}");
  my $backup   = $/;
  my $foldname = $bftr->strip_file_name($file);
  my $sfile    = $bftr->path($settings{'basepath'}, $foldname, "$foldname\.seq");
  my $path     = $bftr->path($settings{'basepath'}, $foldname, "$foldname\.xml");
  my $wrn      = $bftr->path($settings{'basepath'}, $foldname, "warnings\.txt");
  $hash->{'fpath'} = $bftr->path($settings{'basepath'}, "$foldname\/");
  $hash->{'gfile'} = "$path";                             
  $hash->{'sfile'} = "$sfile";
  $hash->{'wfile'} = "$wrn";
  my $xml = "";
  my %w   = ();

  if ( -e $path ) {
    $/ = undef;
    open( IN, $path );
    $xml = <IN>;
    close IN;
    $/ = $backup;
    my $temp = get_xml($xml);
    numerate(\%{$temp->{'gene'}}, 'exon');
    $hash->{'features'} = \%{ $temp->{'gene'} };
    foreach my $exon ( keys %{ $hash->{'features'}{'exon'} } ) {
      my $temp = $hash->{'features'}{'exon'}{$exon};
      numerate($temp, 'warnings');
      $temp->{'from'} = $temp->{'template'}{'from'};
      $temp->{'to'}   = $temp->{'template'}{'to'};
      if ( $temp->{'to'} < $temp->{'from'} ){
        $temp->{'strand'} = -1;
      }
      else {
        $temp->{'strand'} = 1;
      }
      my $c = $temp->{'chromosome'};
      my $f = $temp->{'from'};
      my $t = $temp->{'to'};
      $temp->{'seq'} = $db{"$c:$f,$t"};
    }
  }
  else {
    die("Could not open $path. Please, run bsniffer on this gene\n");
  }
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
        my $n = $bftr->getnkeys( \%{ $hash->{$tag} } );
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
        my $n = $bftr->getnkeys( \%{ $hash->{$tag} } );
        $hash->{$tag}{ $n } = $content;
      }
      else{
        $hash->{$tag} = $content;
      }
    }
  }
  return $hash;
}

sub set_xml {

  # Takes a hash and returns an xml-formatted string with its contents
  # @_ = (hash_ref)
  # This procedure complements get_xml()

  my $hash_ref = shift;
  my $tab      = shift;
  my $prev_key = shift;
  $tab = 1 unless ($tab);
  my $result = '';
  my @keys   = keys %$hash_ref;
  foreach my $key ( sort { lc($a) cmp lc($b) } @keys ) {
    if (ref($hash_ref->{$key}) eq 'HASH') {
      my $text = $hash_ref->{$key};
      if ($key eq '0' || $key =~ /^\d+$/){
        $result .=
          loadTag("$prev_key", set_xml( \%{ $hash_ref->{$key} }, $tab, $prev_key ));
      }
      else{
        my @keys   = keys %{$hash_ref->{$key}};
        my $k = $keys[0];
        if ($k eq '0' || $k =~ /^\d+$/){
          $result .= set_xml( \%{ $hash_ref->{$key} }, $tab, $key )
        }
        else{
          $result .=
            loadTag( "$key", set_xml( \%{ $hash_ref->{$key} }, $tab + 1, $key ));
        }
      }
    }
    else {
      if ($key eq '0' || $key =~ /^\d+$/){
        my $text = $hash_ref->{$key};
        $result .=
          loadTag("$prev_key", $text) if ($text);
      }
      else{
        my $text = $hash_ref->{$key};
        $result .= loadTag( "$key", $text ) if ($text);
      }
    }
  }
  return $result;
}

#sub trim {

#  # Deletes leading and trailing characters in a string
#  # @_ = (string_ref, characters_to_remove);
#  # Returns nothing. Modifies the passed string.

#  my $string = shift;
#  my $todel  = shift;
#  $$string =~ s/^[$todel]+//g;
#  $$string =~ s/[$todel]+$//g;
#}

sub getfromto {

  # Preliminary calculation of the boundaries of the chromosomic
  # region of interest
  # @_ = (hash_ref)
  # Returns nothing. Creates $hash->{'contig'}{'chr'},
  # $hash->{'contig'}{'from'} and $hash->{'contig'}{'to'}
  # This functions can work on unsorted exons. The boundaries
  # include all of the exons plus <delta> nucleotides. Requires
  # a further procedure to make sure that those boundaries are
  # legal (checked in initVS).

  my $hash    = shift;
  my $delta   = shift;
  my $delta_up   = $delta * $hash->{'up_overhang'};
  my $delta_down = $delta * $hash->{'down_overhang'};
  my $mn      = 0;
  my $c       = "";
  my @numbers = ();
  foreach my $n ( keys %{ $hash->{'features'}{'exon'} } ) {
    if ( $hash->{'features'}{'exon'}{$n}{'from'}
      && $hash->{'features'}{'exon'}{$n}{'to'} )
    {
      my $f     = $hash->{'features'}{'exon'}{$n}{'from'};
      my $t     = $hash->{'features'}{'exon'}{$n}{'to'};
      my $frame = ( $f < $t ) ? 1 : -1;
      $f *= $frame;
      $t *= $frame;
      push @numbers, $f;
      push @numbers, $t;
    }
    $c = $hash->{'features'}{'exon'}{$n}{'chromosome'} unless ($c);
  }
  my $f = min(@numbers);
  my $t = max(@numbers);
  $f = abs($f);
  $t = abs($t);
  my $frame = ( $f < $t ) ? 1 : -1;
  if ( $frame == -1 ) {
    $f += $delta_up;
    $t -= $delta_down;
  }
  else{
    $f -= $delta_up;
    $t += $delta_down;
  }
  $c =~ s/>//;
  $hash->{'contig'}{'chr'}  = $c;
  $hash->{'contig'}{'from'} = $f;
  $hash->{'contig'}{'to'}   = $t;
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


sub load_matrix {

  # Loads a probability matrix for Blast operations
  # @_ = (file)
  # Returns a hash containing the assigned probability for
  # every pair in the input matrix
  # At this point, the file must be tab-delimited, with the
  # data in a table format

  my $file   = shift;
  my %result = ();
  open( IN, $file );
  my $line = '#';
  while ( $line && $line eq '#' ) {
    $line = <IN>;
    $line = '#' if ( $bftr->left( $line, 1 ) eq '#' );
  }
  my @hor = split( /\s+/, $line );
  shift @hor;
  while (<IN>) {
    my @fields = split(/\s+/);
    my $fleter = shift @fields;
    for ( my $i = 0 ; $i < @hor ; $i++ ) {
      my $key = uc("$hor[$i]$fleter");
      my $val = $fields[$i];
      $key =~ s/[^\w\d\-\.]//g;
      $val =~ s/[^\w\d\-\.]//g;
      $result{$key} = $val;
    }
  }
  close IN;
  return %result;
}

sub my_blast {

  my $fseq   = shift;
  my $sseq   = shift;
  my $matrix = shift;
  my $result = 0;
  for ( my $i = 0 ; $i < length($fseq) ; $i++ ) {
    my $fl = uc( substr( $fseq, $i, 1 ) );
    my $sl = uc( substr( $sseq, $i, 1 ) );
    if ( exists( $matrix->{"$fl$sl"} ) ) {
      $result += $matrix->{"$fl$sl"};
    }
    else {
      $result -= 0;
    }
  }
  return $result;
}

sub loadBlast {

  ######## Warning: This version of loadBlast only gets hits in the + frames  ###########
  ######## Warning: Template gaps are deleted                                 ###########

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
#  In this version of the routine, the hit '0' would not exist

  my $btext   = shift;
  my $header  = '';
  my $foot    = '';
  my %result  = ();
  my @contigs = $bftr->isplit( '>', $btext );
  $header = shift @contigs unless ( substr( $contigs[0], 0, 1 ) eq '>' );
  for ( my $i = 0 ; $i < scalar @contigs ; $i++ ) {
    my $contig_text = $contigs[$i];
    $contig_text =~ />(.+)\n/;
    my $cname         = $1;
    my %temp_hit_data = ();
    my %hit_data      = ();
    my @hits          = $bftr->isplit( 'Score', $contig_text );
    foreach my $hit (@hits) {
      if ( lc( substr( $hit, 0, length('Score') ) ) eq 'score' ) {
        my %temp_hit = get_hit_data($hit);
        $bftr->hpush( \%temp_hit_data, \%temp_hit ) if ( $bftr->getnkeys( \%temp_hit ) );
      }
    }

    # Sort hits
    foreach my $n (
      sort { $temp_hit_data{$a}{'from'} <=> $temp_hit_data{$b}{'from'} }
      keys %temp_hit_data
      )
    {
      $bftr->hpush( \%{ $hit_data{'hits'} }, \%{ $temp_hit_data{$n} } );
    }
    $bftr->hpush( \%result, \%hit_data );
    $bftr->hadd( \%result, 'contig_name', $cname );
  }
  return %result;
}

sub get_hit_data {

  # Describes a blast hit (see loadBlast)
  # @_ = (hit_text)
  # Returns a hash containing data from the hit_text
  # See loadBlast for an example

  my $text       = shift;
  my %result     = ();
  my $getn       = '.*?([-\+]?\d?e?-?[\d\.]+)';
  my $get_score  = ".*?Score$getn";
  my $get_expect = ".*?Expect$getn";
  my $get_ids    = ".*?Ident$getn$getn$getn";
  $text =~ /$get_score$get_expect$get_ids/si;
  my $score  = $1;
  my $expect = $2;
  my $ids    = $5;
  $ids = $3 / $4 if ( $3 && $4 );

  my $get_seq = '\s([^0-9]+)\s';
  my @lines   =
    $bftr->isplit( "Query$getn$get_seq$getn.*?Sbjct$getn$get_seq$getn", $text );
  my $nl    = @lines;
  my $qtext = '';
  my $ttext = '';
  my $qfrom = 0;
  my $tfrom = 0;
  my $qto   = 0;
  my $tto   = 0;
  if ($nl) {
    for ( my $i = 1 ; $i < $nl ; $i++ ) {
      my $line = $lines[$i];
      $line =~ /Query$getn$get_seq$getn.*?Sbjct$getn$get_seq$getn/si;
      $qfrom = $1 if ( $i == 1 );
      $qto   = $3 if ( $i == $nl - 1 );
      $qtext .= $2;
      $tfrom = $4 if ( $i == 1 );
      $tto   = $6 if ( $i == $nl - 1 );
      $ttext .= $5;
    }
    $qtext =~ s/[^-*\w]//g;
    $ttext =~ s/[^-*\w]//g;
    $result{'score'}            = $score;
    $result{'expect'}           = $expect;
    $result{'ids'}              = $ids;
    $result{'from'}             = $tfrom - 1;
    $result{'to'}               = $tto - 1;
    $result{'template'}{'from'} = $tfrom - 1;
    $result{'template'}{'to'}   = $tto - 1;
    $result{'template'}{'seq'}  = $ttext;
    $result{'query'}{'from'}    = $qfrom;
    $result{'query'}{'to'}      = $qto;
    $result{'query'}{'seq'}     = $qtext;
    $result{'qtext'}            = delinsert( $qtext, $ttext );
  }
  %result = ()
    if ( $result{'to'} < $result{'from'} );    # Empty if the frame is negative
  return %result;
}


sub translate {
  my $seq   = shift;
  my %cdict = (
    'TTT' => 'F',
    'TTC' => 'F',
    'TTA' => 'L',
    'TTG' => 'L',
    'TCT' => 'S',
    'TCC' => 'S',
    'TCA' => 'S',
    'TCG' => 'S',
    'TAT' => 'Y',
    'TAC' => 'Y',
    'TAA' => '*',
    'TAG' => '*',
    'TGT' => 'C',
    'TGC' => 'C',
    'TGA' => '*',
    'TGG' => 'W',
    'CTT' => 'L',
    'CTC' => 'L',
    'CTA' => 'L',
    'CTG' => 'L',
    'CCT' => 'P',
    'CCC' => 'P',
    'CCA' => 'P',
    'CCG' => 'P',
    'CAT' => 'H',
    'CAC' => 'H',
    'CAA' => 'Q',
    'CAG' => 'Q',
    'CGT' => 'R',
    'CGC' => 'R',
    'CGA' => 'R',
    'CGG' => 'R',
    'ATT' => 'I',
    'ATC' => 'I',
    'ATA' => 'I',
    'ATG' => 'M',
    'ACT' => 'T',
    'ACC' => 'T',
    'ACA' => 'T',
    'ACG' => 'T',
    'AAT' => 'N',
    'AAC' => 'N',
    'AAA' => 'K',
    'AAG' => 'K',
    'AGT' => 'S',
    'AGC' => 'S',
    'AGA' => 'R',
    'AGG' => 'R',
    'GTT' => 'V',
    'GTC' => 'V',
    'GTA' => 'V',
    'GTG' => 'V',
    'GCT' => 'A',
    'GCC' => 'A',
    'GCA' => 'A',
    'GCG' => 'A',
    'GAT' => 'D',
    'GAC' => 'D',
    'GAA' => 'E',
    'GAG' => 'E',
    'GGT' => 'G',
    'GGC' => 'G',
    'GGA' => 'G',
    'GGG' => 'G',
    'TCN' => 'S',
    'CGN' => 'R',
    'CCN' => 'P',
    'CTN' => 'L',
    'ACN' => 'T',
    'GGN' => 'G',
    'GCN' => 'A',
    'GTN' => 'V',
  );
  my $result = '';
  my @codons = unpack( '(A3)*', $seq );
  foreach my $codon (@codons) {
    if ( exists( $cdict{uc($codon)} ) ) {
      $result .= $cdict{uc($codon)};
    }
    else {
      $result .= "X";
    }
  }
  return $result;
}


sub pformat {

  # Formats a protein sequence for display
  # @_ = (prot_seq)
  # Returns a string with spaces for a protein sequence to fit its DNA sequence.

  my @chars = split( //, $_[0] );
  my $qseq = "";
  foreach my $c (@chars) {
    $qseq .= "$c  ";
  }
  return $qseq;
}

sub ungap {

# Converts a relative position in a gapped sequence into a relative position in the
# corresponding ungapped sequence (i. e. the same position if we do not count the gaps).
# @_ = (seq, pos)
# <seq> -> gapped sequence, string.
# <pos> -> query position in <seq>, counting gaps, integer.
# Returns an integer corresponding to <pos> if gaps are not included in <seq>.

  my $seq  = shift;
  my $upos = shift;
  return $upos if ( $upos <= 0 );
  my @a = split( //, $seq );
  my $result = 0;
  for ( my $i = 0 ; $i < $upos ; $i++ ) {
    if ( !$a[$i] || $a[$i] ne '-' ) {
      $result++;
    }
  }
  return $result;
}

sub gap {

# Converts a relative position in an ungapped sequence into a relative position in the
# corresponding gapped sequence (i. e. the same position if we count the gaps).
# @_ = (seq, pos)
# <seq> -> gapped sequence, string.
# <pos> -> query position in <seq>, not counting gaps, integer.
# Returns an integer corresponding to <pos> if gaps are included in <seq>.

  my $seq  = shift;
  my $upos = shift;
  return $upos if ( $upos <= 0 );
  my @a = split( //, $seq );
  my $result = 0;
  for ( my $i = 0 ; $i < $upos ; $i++ ) {
    if ( $a[$i] eq '-' ) {
      $result++;
    }
    $result++;
  }
  return $result;
}

sub loadTag {
  my ( $tid, $ttext ) = @_;
  $ttext = '1' unless ($ttext);
  chomp $ttext;
  $ttext =~ s/(\n)/\n\t/gs;
  my $result = "<$tid>\n\t$ttext\n<\/$tid>\n";
  return $result;
}


#########################Debug######################
sub hash_to_string {
  my $hash_ref = shift;
  my $tab      = shift;
  $tab = 1 unless ($tab);
  my $result = '';
  my @keys   = keys %$hash_ref;
  foreach my $key ( sort { lc($a) cmp lc($b) } @keys ) {
    if ( ref( $hash_ref->{$key} ) eq 'HASH' ) {
      my $text = $hash_ref->{$key};
      my $t = sp($tab) . sp( length("$key => "), " " );
      $text =~ s/\n/\n$t/g;
      $result .=
          sp($tab)
        . "$key => $text\n"
        . hash_to_string( \%{ $hash_ref->{$key} }, $tab + 1 );
    }
    elsif ( ref( $hash_ref->{$key} ) eq 'ARRAY' ) {
      $result .=
        sp($tab) . "$key =>" . join( ', ', @{ $hash_ref->{$key} } ) . "\n";
    }
    else {
      my $text = $hash_ref->{$key};
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

package GT_Term_Screen;
########################Screen management###########

eval 'require Term::Screen';

sub new{
  my $pkg = shift;
  my $self = {};
  bless($self, $pkg);
  return $self;
}

sub tprompt {

  # Displays a message and waits for a one-key answer from the user
  # @_ = (master_hash, prompt, possibilities)
  # <master_hash>   -> hash with all the information, hash reference
  # <prompt>        -> text to display, string
  # <possibilities> -> possible answers, string. By default, 'yn'.
  # Returns a string with the answer of the user, only if it is one
  # of the characters in <possibilities>, case independent.

  my $self = shift;
  my $hash = shift;
  my $msg  = shift;
  my $poss = shift;
  $poss = 'yn' unless ($poss);
  $hash->{'warnings'} = {};
  $self->message($hash, $msg);
  $self->showVS($hash);
  my $answer = '';
  $poss = "\[$poss\]";

  while ( !( $answer =~ /$poss/i ) ) {
    $answer = lc( $self->get_char() );
  }
  return $answer;
}

sub prompt_for_expression {

  # Warning: this function contains interface-specific commands
  # @_ = (master_hash, prompt)
  # <master_hash>   -> hash with all the information, hash reference
  # <prompt>        -> text to display, string
  # Returns a string with the answer of the user

  my $self = shift;
  my $hash   = shift;
  my $disp   = shift;
  my $result = "";
  $hash->{'warnings'} = {};
  $self->message($hash, $disp);
  $self->showVS($hash);
  $self->endScr();
  $result = vBIO::BATI->prompt();
  $self->initScr();
  return $result;
}

sub initScr {

  # Initializes the interface. Can be extended for several interfaces with
  # a global setup variable. At this point, it is specific for Term::Screen.

  my $self = shift;
  if ( !exists($self->{'scr'}) || !$self->{'scr'} ) {
    $self->{'scr'} = eval 'new Term::Screen';
    $self->{'scr'}->clrscr();
    $self->{'scr'}->flush_input();
    $self->{'scr'}->at( 0, 0 );
  }
}

sub endScr {

  # Ends the interface.

  my $self = shift;
  delete($self->{'scr'}) if (exists($self->{'scr'}));
}

sub get_char {
  my $self = shift;
  my $result = $self->{'scr'}->getch();
  return $result;
}

sub showVS {
  my $self = shift;
  my $hash = shift;
  $self->{'scr'}->clrscr();
  $self->{'scr'}->resize();
  my $nrows = $self->{'scr'}->rows();
  my $ncols = $self->{'scr'}->cols();
  main::refreshS( $hash, $ncols, $nrows );
  foreach my $n ( keys %{ $hash->{'screen'} } ) {
    $self->{'scr'}->at( $n, 0 );
    $self->{'scr'}->puts( $hash->{'screen'}{$n} );
  }
  my %hl = main::decodehl( \%{ $hash->{'highlights'} } );
  foreach my $n ( keys %hl ) {
    if ( $hl{$n}{'effect'} == 2 ) {
      $self->{'scr'}->bold();
      $self->{'scr'}->reverse();
    }
    elsif ( $hl{$n}{'effect'} == 1 ) {
      $self->{'scr'}->bold();
    }
    else {
      $self->{'scr'}->reverse();
    }
    $self->{'scr'}->at( $hl{$n}{'y'}, $hl{$n}{'x'} );
    $self->{'scr'}->puts( $hl{$n}{'text'} );
    $self->{'scr'}->normal();
  }
  $hash->{'highlights'} = ();
  $self->{'scr'}->at( $nrows - 1, 0 );
}

sub message{
  my $self = shift;
  my $hash = shift;
  my $msg  = shift;
  vBIO::BATI->hpush(\%{$hash->{'warnings'}}, $msg);
}

1;

package GT_Win32_Screen;

########################Screen management###########

eval 'use Win32::Console';


sub new{
  my $pkg = shift;
  my $self = {};
  bless($self, $pkg);
  return $self;
}

sub tprompt {

  # Displays a message and waits for a one-key answer from the user
  # @_ = (master_hash, prompt, possibilities)
  # <master_hash>   -> hash with all the information, hash reference
  # <prompt>        -> text to display, string
  # <possibilities> -> possible answers, string. By default, 'yn'.
  # Returns a string with the answer of the user, only if it is one
  # of the characters in <possibilities>, case independent.

  my $self = shift;
  my $hash = shift;
  my $msg  = shift;
  my $poss = shift;
  $poss = 'yn' unless ($poss);
  $hash->{'warnings'} = {};
  $self->message($hash, $msg);
  $self->showVS($hash);
  my $answer = '';
  $poss = "\[$poss\]";

  while ( !( $answer =~ /$poss/i ) ) {
    $answer = lc( $self->get_char());
  }
  return $answer;
}

sub prompt_for_expression {

  # Warning: this function contains interface-specific commands
  # @_ = (master_hash, prompt)
  # <master_hash>   -> hash with all the information, hash reference
  # <prompt>        -> text to display, string
  # Returns a string with the answer of the user

  my $self = shift;
  my $hash   = shift;
  my $disp   = shift;
  my $result = "";
  $hash->{'warnings'} = {};
  $self->message($hash, $disp);
  $self->showVS($hash);
#  endScr();
  $result = vBIO::BATI->prompt();
  $self->initScr();
  return $result;
}

sub message{
  my $self = shift;
  my $hash = shift;
  my $msg  = shift;
  vBIO::BATI->hpush(\%{$hash->{'warnings'}}, $msg);
}


sub initScr {

  # Initializes the interface. Can be extended for several interfaces with
  # a global setup variable.

  my $self = shift;

  if ( !exists($self->{'scr_in'}) || !$self->{'scr_in'} ) {
    $self->{'scr_in'} = eval'new Win32::Console(Win32::Console::STD_INPUT_HANDLE)';
    if (!$self->{'scr_in'}){
      die("Problem with the input interface: $@\n");
    }
    #$self->{'scr_in'}->Cls();
    #$self->{'scr'}->flush_input();
    #$self->{'scr_in'}->at( 0, 0 );
  }
  if ( !$self->{'scr_out'} ) {
    $self->{'scr_out'} = eval'new Win32::Console(Win32::Console::STD_OUTPUT_HANDLE)';
    if (!$self->{'scr_out'}){
      die("Problem with the output interface: $@\n");
    }
    $self->{'scr_out'}->Cls();
    $self->{'scr_out'}->Cursor( 0, 0 );
  }
}

sub endScr {

  # Ends the interface. Not necessary with Win32::Console.
  # Kept for consistency.

#  $self->{'scr_in'}  = "" if ($self->{'scr_in'});
#  $self->{'scr_out'} = "" if ($self->{'scr_out'});
}

sub get_char{
  my $self = shift;
  my $result = '';
  my @event = $self->{'scr_in'}->Input();
  if($event[0] && $event[0] == 1 && $event[1]) {
      # LEFT ARROW
      if( $event[3] == 37
        && $event[4] == 75){
        $result = 'kl';
      }
      # RIGHT ARROW
      elsif( $event[3] == 39
        && $event[4] == 77){
        $result = 'kr';
      }
      # UP ARROW
      elsif( $event[3] == 38
        && $event[4] == 72){
        $result = 'ku';
      }
      # DOWN ARROW
      elsif( $event[3] == 40
        && $event[4] == 80){
        $result = 'kd';
      }
      # F1
      elsif( $event[3] == 112){
        $result = 'k1';
      }
      # F2
      elsif( $event[3] == 113){
        $result = 'k2';
      }
      # F3
      elsif( $event[3] == 114){
        $result = 'k3';
      }
      else{
        $result = chr($event[5]);
      }
  }
  return $result;
}


sub showVS {
  my $self = shift;
  my $hash = shift;
  my $cursor = 0;
  #$self->{'scr'}->resize();
  $self->{'scr_out'}->Cls();
  my ($wLeft, $wTop, $wRight, $wBottom) = $self->{'scr_out'}->Window();
  my $nrows = $wBottom - $wTop;
  my $ncols = $wRight - $wLeft;
  main::refreshS( $hash, $ncols, $nrows );
  foreach my $n ( keys %{ $hash->{'screen'} } ) {
    $self->{'scr_out'}->Cursor( 0, $n );
    $self->{'scr_out'}->Write( $hash->{'screen'}{$n} );
  }
  my %hl = main::decodehl( \%{ $hash->{'highlights'} } );
  foreach my $n ( keys %hl ) {
    my $l = length($hl{$n}{'text'});
    my $attr_bold = eval 'chr($Win32::Console::FG_RED) x $l';
    my $attr_rev  = eval 'chr($Win32::Console::ATTR_INVERSE) x $l';
    my $attr_bold_rev  = eval 'chr($Win32::Console::FG_RED|$Win32::Console::ATTR_INVERSE) x $l';
    $self->{'scr_out'}->Cursor( $hl{$n}{'x'}, $hl{$n}{'y'} );
    $self->{'scr_out'}->Write( $hl{$n}{'text'} );
    if ( $hl{$n}{'effect'} == 2 ) {
      $self->{'scr_out'}->WriteAttr($attr_bold_rev, $hl{$n}{'x'}, $hl{$n}{'y'});
    }
    elsif ( $hl{$n}{'effect'} == 1 ) {
      $self->{'scr_out'}->WriteAttr($attr_bold, $hl{$n}{'x'}, $hl{$n}{'y'});
    }
    else {
      $self->{'scr_out'}->WriteAttr($attr_rev, $hl{$n}{'x'}, $hl{$n}{'y'});
    }
  }
  $hash->{'highlights'} = ();
  $self->{'scr_out'}->Cursor( 0, $nrows - 1 );
}


1;



=head1 NAME

GeneTuner, a script for the manual curation of gene
families.

=head1 SYNOPSIS

genetuner [project_name] [gene_name]

=head1 DESCRIPTION

This program allows the user to define the exons of a set of
genes in a set of genomic contigs. It is particularly suited
for similarity-based predictions applied to gene families.
Usually, GeneTuner is run downstream of BlastSniffer, which
provides all of the input files this program needs. The
purpose of both programs is to automatically perform several
time-consuming steps in the manual annotation of genes by
homology.

Usually, annotations start with a I<TBLASTN> comparison
between a known protein and a genomic sequence. The user
must then decide which resulting hits correspond to I<bona
fide> exons of new putative genes. BlastSniffer handles this
task. Then, every exon junction has to be manually assessed
and corrected, and the new gene must be rebuilt by copying
the specified subsequences from the template contig.
GeneTuner lets the user choose exon junctions, and also add
and remove exons, directly from the genomic sequence. To
help this task, the template sequence is shown aligned to
the contig.

Furthermore, the user does not need to keep track of the
frame. The length of every exon is a multiple of three, so
that the reading frame is always correct. At the end, the
program decides automatically if each junction must be
edited to improve the similarity to the template sequence
and/or to get canonical exon/intron junctions.

The recommended procedure to use these programs consists of
the following steps:

=over

=item 1. Store all of the starting protein sequences in one
folder (F<ppath>) with a given extension (F<aa_ext>).
Optionally, store the corresponding nucleotide sequences
with the same names in a different folder (F<npath>) and/or
with a different extension (F<nt_ext>).

=item 2. Run F<tbex> and create a new project (see tbex
documentation). This will compare all of the sequences in F<ppath>
with the target genomic sequence by tblastn. Optionally, at
this step you can also run F<bgmix> on the same project. This
will create a composite file with all of the tblastn results
which may be useful to call difficult genes, like tandem repeats.

=item 3. Run BlastSniffer on those files (see bsniffer
documentation).

=item 4. Open GeneTuner and load the project file created by
BlastSniffer. This can be done from the command line or by
choosing I<Load a project> at the welcome screen. Every time
you load a project, a different F<tbn> file is edited,
unless I<gene_name> is provided after I<project_name>.

=item 5. Edit the predicted exons, save the result, and
exit. Repeat steps 4 and 5 until you have processed every
F<tbn> file.

=item 6. Run F<bgmix> on the project. This will create a composite
file with all of the tblastn comparisons and highlight those hits
overlaping defined exons. This helps the identification of
novel putative genes which have not been annotated. To annotate
these genes, you can add suitable protein sequences to F<ppath>
and start at step 1. You can also add the novel proteins from
the previous steps to search for paralogues.

=back

=head2 Welcome screen

Unless a valid project name is provided in the command line,
the user will be taken to the welcome screen. This screen
allows the user to create a new project (usually not
necessary) or load an existing project. Once a project is
created or loaded, GeneTuner processes the next F<tbn> file
and takes the user to the edition screen.

=head2 Edition screen

These are the elements of the edition screen:

  ____________________________________________________
 |                                                    |
 |c19h_usp1       length: 785                         |(1)
 |Chromosome: Chr8        strand: 1                   |(2)
 |                                                    |
 |                                                    |
 |                                                    |
 |                                                    |
 |                  |25424819                         |(3)
 |                           995            1009      |(4)
 |                           gttcctgcagcacag          |(4)
 |                           |||||||||||||||          |
 | TTTCTCTGTTTCCCAGTGACCAGGTGGTTCCTGCAGCACAGCCCCCCTC  |(5)
 |                                                    |
 | F  L  C  F  P  V  T  R  W  F  L  Q  H  S  P  P  P  |(6)
 |  F  S  V  S  Q  *  P  G  G  S  C  S  T  A  P  L    |(6)
 |   S  L  F  P  S  D  Q  V  V  P  A  A  Q  P  P  S   |(6)
 |                  |  |  |  |  |  |  |  |        |   |
 |                  D  Q  V  V  P  A  A  Q  S  -  S   |(7)
 |                  58                                |(7)
 |                                                    |
 |                                                    |
 |                                                    |
 |                                                    |
 |Lacking residues                                    |(8)
 |____________________________________________________|

=over

=item 1. Name of the F<tbn> file in process, and length of
the protein sequence.

=item 2. Name and strand (+/-) of the chromosome (or contig)
beeing displayed.

=item 3. Current cursor position in the contig.

=item 4. Matching sequence in the starting nucleotide
sequence, found by I<BLASTN>.

=item 5. Sequence of the contig.

=item 6. Translation of the contig in all three frames in
the displayed I<strand>.

=item 7. Matching sequence in the starting protein sequence,
found by I<TBLASTN>.

=item 8. Warnings related to the current exon. Messages and
prompts for the user are also displayed here.

=back

=head2 Controls

The editing screen is controlled with the keyboard. Controls
are case-insensitive.

=over

=item * Left and right arrows scroll the contig one base at
a time. If the user is marking an exon, the contig is
scrolled three bases at a time.

=item * "O" and "P" scroll the contig a number of bases at a
time. This number is 1000 by default, but can be configured.

=item * "M" gets the cursor to the previous "atg" sequence.

=item * Space bar starts or ends an exon at the current
cursor point.

=item * "Z" / "X" extend the previous / next exon up to the
current cursor point.

=item * "A" deletes the exon the cursor is located in.

=item * "W" adds a warning to the current exon. If that
warning is left empty, deletes current warnings.

=item * Up and down arrows get the cursor to the
previous/next start or end of an exon.

=item * "K" / "L" get the cursor to the previous / next
start or end of a I<BLASTN> hit.

=item * "," / "." get the cursor to the previous / next
start or end of a I<TBLASTN> hit.

=item * "R" / "F" increase / decrease the limit expectation
value (E) of the internal BLASTN search (4).

=item * "E" / "D" increase / decrease the limit expectation
value (E) of the internal TBLASTN search (7).

=item * "S" performs a local TBLASTN search with automatic
parameters.

=item * "(" / ")" sets boundaries for local TBLASTN search.

=item * "Y" / "H" increase / decrease the limit expectation
value (E) of the current local TBLASTN search.

=item * "Ctrl+F" sets a nucleotide pattern to search. See
I<Advanced search>.

=item * "Ctrl+B" sets a peptide pattern to search. See
I<Advanced search>.

=item * F2 finds current pattern backward in the contig.

=item * F3 finds current pattern forward in the contig.

=item * "Ctrl+G" takes the cursor to a given position in the
contig.

=item * "Ctrl+Z" undoes the last exon modification.

=item * "E<lt>" reads an additional 5000 bases chunk
upstream from the contig shown.

=item * "E<gt>" reads an additional 5000 bases chunk
downstream from the contig shown.

=back

=head2 Navigation

Basic navigation through the contig is performed with the
cursor keys. Left and right move the cursor one base at a
time in normal mode and three bases at a time in exon
extension mode. Up and down arrows place the cursor in the
previous or next exon boundary. Since the alignment of the
contig with the template sequence is very important, it can
be browsed like the exons, with "K" and "L" for the
I<BLASTN> alignment and "," and "." for the I<TBLASTN>
alignment.

=head2 Edition of exons

When a putative gene is loaded, its genomic sequence is
shown, with 5000 extra bases upstream and downstream. If
more extra bases are needed, E<lt> adds 5000 extra bases
upstream and E<gt> adds 5000 extra bases downstream.

The predicted exons are highlighted. To start a new exon,
place the cursor at the beginning and press the space bar.
Then, extend the exon up to the end and press the space bar
again. While extending the exon, the minimal step is changed
to three bases instead of one, so that the translation frame
is always conserved.

If the space bar is pressed while the cursor is inside an
exon, the user will enter extension mode for that exon. To
directly change the start or end of an exon, place the
cursor on the new boundary and press "Z" to change the start
of the next exon or "X" to change the end of the previous
exon. If necessary, the boundary will be automatically
shifted to preserve the translation frame.

To remove an exon, simply place the cursor inside that exon
and press "A". The last exon edition can be undone by
pressing Ctrl+Z.

=head2 Advanced search

When looking for dissimilar and/or short exons, it may be
necessary to use some advanced search capabilities. First,
the limit expect value of the I<BLASTN> and I<TBLASTN>
comparisons may be lowered to improve their sensitivity.
However, this will also increase the noise in the BLAST
comparison. To avoid this noise increase, local I<TBLASTN>
comparisons can be performed with selected chunks of query
protein and genomic sequence. Finally, users can perform
nucleotide and peptide pattern searches.

=head3 BLAST sensitivity

By default, GeneTuner performs gapped, unfiltered BLASTN and
TBLASTn comparisons with a limit expect value of 10. This
value can be raised or lowered with the keys "E" and "D".
Raising the limit expect value is useful when long stretches
of the starting and target sequences are dissimilar. On the
other hand, this procedure is likely to produce multiple
non-significant hits, which may also slow down execution of
the program.

=head3 Local TBLASTN

TBLASTN searches can be confined to specific regions of the
query protein and the target contig. This prodecure allows
the user to increase the sensitivity of the search while
avoiding hits outside the region of interest. Also, this
mode is optimized for short query sequences. The boundaries
for this search can be defined in two complementary ways:

=over

=item * Pressing "S" prompts the program to suggest hits to
fill a void in the annotation. The program locates the
previous and next exons from the cursor, which are set as
genomic boundaries. The corresponding query protein
boundaries are calculated based on the internal TBLASTN
comparison. Then, a local TBLASTN comparison is performed
with the protein sequence and the genomic sequence between
these calculated boundaries. If no new hits are found, the
comparison is repeated with increasing limit expect values.
This method is fast, but less reliable than setting the
boundaries manually.

=item * Boundaries can also be set manually. "(" sets the
start genomic boundary at the current cursor position, while
")" sets the finish genomic boundary. The user is prompted
to enter the corresponding query protein boundaries.

=back

The limit expect value for the current local TBLASTN
comparison can be raised or lowered by pressing "Y" and "H".
Local hits are discarded when another local or general
TBLASTN is performed.

=head3 Search patterns

The genomic contig can be probed with nucleotide and
peptidic patterns. This allows the user to find very short
exons and conserved motifs which may not be amenable to
BLAST searches.

To enter a nucleotide pattern, press Ctrl-F. To enter a
peptide pattern, press Ctrl-B. Then, press F2 to find the
previous occurrence of the pattern upstream from the cursor
or F3 to find the next occurrence of the pattern downstream
from the cursor. Patterns are based on regular expressions
from Perl. The following elements are supported:

=over

=item * "." means "any nucleotide". "X" means "any residue".

=item * "[ACD]" means "either A, C or D" nucleotide or
residue.

=item * "*" means "0 or more occurrences".

=item * "{n}" means "exactly n occurrences". {n,} means "n
or more occurrences". {n, m} means "n or more occurrences up
to m".

=item * "Z" means a stop codon (only in peptide mode).

=back



=head2 Project files

Project files are saved into the F<projects> folder with
F<.gt> extension. They contain the following information:

 tbnpath  => path to the folder where the starting ".tbn"
             files are stored
 basepath => path to the folder where the result folders
             will be stored
 dbpath   => path to the genomic (template) sequence
 ppath    => path to the starting protein sequences folder
 aa_ext   => extension of the starting protein sequences
 npath    => path to the starting nucleotide sequences
             folder (optional)
 nt_ext   => extension of the starting nucleotide sequences
             (optional)

Both GeneTuner and BlastSniffer can create compatible
project files from user input.

=head1 ARGUMENTS

 --help        print this help

 =head1 OPTIONS

none implemented

=head1 FILES

 genetuner.pl       script file

=head1 DEPENDENCIES

GeneTuner requires BioPerl, namely 'Bio::DB::Fasta' and
'Bio::SeqIO' modules. Win32 version also requires
Win32::Console, while *NIX version requires Term::Screen.

=head1 INPUT

BlastSniffer creates input files GeneTuner needs. For a gene
called F<example.aa>, GeneTuner needs one file called
F<example.xml> at the F<basepath/example> folder. This file
must contain an extensible markup format (XML) with the
description of the predicted gene. Here is a simple example:

 <gene>
   <exon>
     <chromosome>chr8</chromosome>
     <template>
       <from>28268822</from>
       <to>28268947</to>
     </template>
     <warnings>Query does not start at 1</warnings>
   </exon>
   <exon>
     <chromosome>chr8</chromosome>
     <template>
       <from>28269898</from>
       <to>28270059</to>
     </template>
   </exon>
 </gene>


=head1 OUTPUT

=over

=item * An F<.xml> file with annotation info.

=item * A F<.seq> file with the predicted sequence, one exon
per line.

=item * A file called F<warnings.txt> with a summary of user
notes.

=item * A file called F<map.html> with a genomic map of the
gene. Large introns can be hidden by clicking on them. To
show these introns again, click on the horizontal line that
appears.

=item * A F<.gff> file with the annotation in GFF format.
This file can be used as input to other visualization and
edition programs, like gff2ps or UCSC and Ensembl genomic
browsers.

=back

=head1 LICENSE

This program is free software and can be redistributed under
the same terms as Perl. See
http://www.perl.com/pub/a/language/misc/Artistic.html

=head1 AUTHOR

Copyright (C) 2008, Victor Quesada

e-mail: quesadavictor@uniovi.es

=cut
