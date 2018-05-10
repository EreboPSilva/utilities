#!/usr/bin/perl -w
use strict;
use vDOM::ps;


my $defs = <<END;
/fsize 6 def
/ypadding 40 def
/width 500 def
%newpath
%width 700 moveto
%width 100 lineto
%0 0 0 setrgbcolor
%stroke
/height 700 def
/y 750 def
/x 150 def
/cxpad 1.7 def % Horizontal padding per character
/cypad 2.3 def % Vertical padding per character
/ydeb 100 def
/names_in_every_line true def

%Fonts

/namefont {
  /Arial-Bold findfont                                                                                    
  fsize scalefont                                                                                         
  setfont
} bind def

/seqfont  {
  /Arial findfont                                                                                    
  fsize scalefont                                                                                         
  setfont
} bind def

/minid 0.7 def   % Light blue for positions with a higher conservation
                 % Dark blue for completely conserved positions
/color_background [1 1 1] def
/color_font [0 0 0] def
/color1 [0.2 0.2 1] def
/color2 [0.7 0.7 1] def

/arrow{
	/alpha_R 150 def
	/afill true def
	/y2 exch def
	/x2 exch def
	/y1 exch def
	/x1 exch def
	/modul x2 x1 sub dup mul y2 y1 sub dup mul add sqrt def
	/R modul 2 div def
	/alpha alpha_R R div def
	/cosa alpha cos def
    /sina alpha sin def
	/dy y2 y1 sub def
	/dx x2 x1 sub def
	/l dx dx mul dy dy mul add sqrt def
	/cosb dx l div def
	/sinb dy l div def
	/coscos cosa cosb mul def
	/sinsin sina sinb mul def
	/sincos sina cosb mul def
	/cossin cosa sinb mul def
	%x3%
	x2
	coscos sinsin add
	R mul
	sub
	/x3 exch def
	%y3%
	y2
	sincos cossin sub
	R mul
	add
	/y3 exch def
	%x4%
	x2
	coscos sinsin sub
	R mul
	sub
	/x4 exch def
	%y4%
	y2
	cossin sincos add
	R mul
	sub
	/y4 exch def

	newpath
		x1 y1 moveto
		x2 y2 lineto
		0 0 0 setrgbcolor
		1 setlinewidth
	stroke
	newpath
		x3 y3 moveto
		x2 y2 lineto
		x4 y4 lineto
		0 0 0 setrgbcolor
		1 setlinewidth
		1 setlinejoin
		afill{
			gsave
			closepath
			fill
			grestore
		} if
	stroke
} def

/mtext_size{
  /text exch def
  newpath
  0 0 moveto
  text false charpath 
  pathbbox
  /th exch def
  /tw exch def
  pop pop
  tw th
} def

/debug{
  30 string cvs print
  (\n) print
} def

/pusharr{
  /val exch def
  /paarr exch def
  /n paarr length 1 add def
  /result n array def
  0 1 paarr length 1 sub{
    /i exch def
    result i paarr i get put
    
  } for
  result n 1 sub val put
  result
} def


/debugarr {
  /arr exch def
  arr {
    print
  } forall
  (
) print flush
} def

/textshade {
  /tscol exch def
  /tsy exch def
  /tsx exch def
  newpath
  tsx cw 2 div sub cy ch 4 div sub moveto
  0 ch rlineto
  cw 0 rlineto
  0 ch neg rlineto
  closepath
  tscol aload pop setrgbcolor
  gsave
  fill
  grestore
  0.5 setlinewidth
  stroke
} def


/getfreq{
  /gfres exch def
  /gfpos exch def
  /gfarr exch def
  /gfresult 0 def
  gfres (-) ne{
	  /gftotal gfarr length def
	  0 1 gftotal 1 sub {
		/gfi exch def
		/gftres gfarr gfi get gfpos get def
		gftres gfres eq {
		  /gfresult gfresult 1 add def
		} if
	  } for
	  /gfresult gfresult gftotal div def
  } if
  gfresult
} def

/write_names {
  namefont
  /tmpy cy def
  0 1 nseq 1 sub {
    /k exch def
	/name seqnames k get def
	name mtext_size pop
	/snw exch def
	/myx cx snw sub cw 1 mul sub def
	myx cy moveto name show
	/cy cy ch sub def
  } for
  /cy tmpy def
  seqfont
} bind def

/next_line{
  /cx x def
  /cy cy nseq ch mul sub 2 ch mul sub def
  /nextcy cy ch nseq mul sub def
  nextcy y height sub lt{
    grestore
    showpage
    gsave
    /cy y def
  } if
  names_in_every_line{
    write_names
  } if
} def

/addif{
  /aiai exch def
  aiai (-) eq {
    0
  } {
    1
  } ifelse
} def

/tostring{
  30 string cvs
} def

/highlight_this{
  cx cy ch 2 mul add cx cy ch add arrow
  seqfont
  arpos tostring mtext_size
  /htth exch def
  /httw exch def
  cx httw 2 div sub cy htth 2 div add  ch 2 mul add moveto arpos tostring show
  
} def

/alig_render{
  /to exch def
  /from exch def
  /high exch def
  /arr exch def
  /cy y def
  /cx x def
  /nseq arr length def
  /seqnames 0 array def  % Array with names
  /seqs 0 array def      % Array with arrays of chars (seqs)
  /str1 arr 0 get 1 get def
  /alig_length str1 length def
  /arhln 0 def           % Current highlight array element
  /arhle high arhln get def
  /hln high length def
  % Pre-processing
  arr{
    /el exch def
	/seqname el 0 get def
	/seq el 1 get def
	/seqnames seqnames seqname pusharr def
	/tmp alig_length array def
	0 1 alig_length 1 sub {
	  /i exch def
	  /char seq i 1 getinterval def
	  tmp i char put
	} for
	/seqs seqs tmp pusharr def
  } forall
  % Set seqnames
  namefont
  write_names
  % Set seqs
  seqfont
  0 1 alig_length 1 sub {
    % Should I print this column?
	/printme arpos 2 add from gt arpos to lt and def  
    /i exch def
	  cx width gt{
	    next_line
	  } if
	  arpos arhle eq{   % highlight this one
	    printme {highlight_this} if
	    % Get next highlight?
	    arhln hln 1 sub lt{
	      /arhln arhln 1 add def
	      /arhle high arhln get def
	    } if
	  } if
	  /arpos arpos seqs 0 get i get addif add def
	  printme{
	      /cx cx cw add def
		  /mytmpy cy def
		  0 1 nseq 1 sub {
			/freq 0 def
			/j exch def
			/res seqs j get i get def
			/mycol color_background def
			/freq seqs i res getfreq def
			freq 1 eq {
			  /mycol color1 def
			} {
			  freq minid gt{
				/mycol color2 def
			  } if 
			} ifelse
			cx cy mycol textshade
			/tw res mtext_size pop 2 div def
			color_font aload pop setrgbcolor
			cx tw sub cy moveto res show
			/cy cy ch sub def
		  } for
		  /cy mytmpy def
	  } if
  } for
} def

seqfont
(M) mtext_size
/mth exch def
/mtw exch def
/cw mtw cxpad mul def   % Character height plus padding
/ch mth cypad mul def   % Character width plus padding
% Draw


END


my $extra = shift;   # Sequences in fasta format
my @temp = split /[\/\\]/, $0;
my $pname = pop @temp;
die("Use: perl $pname multifasta_file\n")
  unless ($extra && -e $extra);


my @fastas = ();
my $x = slurp($extra);


my $alig = `clustalo -i $extra`;

push @fastas, $alig;

foreach my $fasta (@fastas){
  my ($name) = $fasta =~ />(\S+)/;
  warn ("$name\n");
  my $ps = vDOM::ps->new(1);
  $ps->addlines($defs);
  
  my $f = getfasta->new();
  $f->add_seqs($fasta);
  my $start = 1;
  my $end = $f->seqlen(0);
  my $n = $f->convert_pos_to_alig($start);
  my $m = $f->convert_pos_to_alig($end);
  my $string = $f->prepare_for_alignment('all', $n, $m);
  
  #my $start = $n + 1;
  $ps->addlines("/arpos $start 1 sub def",
                "/startpos arpos tostring def", 'startpos mtext_size',
                '/sth exch def', '/stw exch def',
                'x cw add stw 2 div sub y ch add sth 1 mul add moveto',
                '%startpos show');
  $ps->addlines('/alig [', $string, '] def');
  my $markers = '';
  for (1 .. 50){
    my $t = 20 * $_;
    $markers .= "$t ";
  }
  $ps->addlines("/highlight [$markers] def");
  $ps->addlines('alig highlight 1 10000 alig_render');
  my $fname = "$extra\.ps";
  print "$fname\n";
  open(OUT, ">$fname");
  print OUT $ps->get_code();
  close OUT;
}


sub slurp{
  my $file = shift;
  my $result = '';
  local $/;
  open (IN, $file) or die ("Could not open $file: $!\n");
  $result = <IN>;
  close IN;
  return $result;
}



package getfasta;

sub new{
  my $pkg   = shift;
  my $ffile = shift;
  my $ext   = shift || '';
  my $self = {'file' => $ffile, 'errors' => []};
  bless($self, $pkg);
  return $self;
}

sub add_file{
  my $self = shift;
  my $file = shift;
  if (!-e $file){
    push @{$self->{'errors'}}, "$file does not exist";
  }
  $self->{'file'} = $file;
  open (IN, $self->{'file'}) or die($self->{'file'});
  local $/;
  my $text = <IN>;
  $self->_init($text);
  return 1;
}

sub add_seqs{
  my $self = shift;
  my $text = shift;
  $self->_init($text);
  return 1;
}

sub add_folder{
  my $self = shift;
  my $folder = shift;
  if (substr($folder, -1, 1) ne '/'){
    $folder .= '/';
  }
  my $ext   = shift || '';
  return 0 unless (-d $folder);
  opendir DIR, $folder;
  while (my $f = readdir(DIR)){
    if ($ext){
      next unless ($f =~ /\.$ext$/);
    }
    $self->add_file($folder.$f);
  }
  closedir DIR;
}

sub _init{
  my $self = shift;
  my $text = shift;
  while ($text =~ /^>([^\s;]+).*?\n([^>]+)/gsm){
    my $name = $1;
    my $seq = $2;
    $seq =~ s/\s//gsm;
    push @{$self->{'seqs'}}, {'name' => $name, 'seq' => uc($seq)};
    $self->{'byname'}{$name} = uc($seq);
  }
  close IN;
  return 0;
}

sub get_codon{
  my $self = shift;
  my $nseq = shift;
  my $pos  = shift;
  if (!defined($pos)){
    $pos = $nseq;
    $nseq = 0;
  }
  my $frame = $pos % 3;
  my $start = $pos - $frame;
  my $result = substr($self->{'seqs'}[$nseq]{'seq'}, $start, 3);
  return $result;
}

sub get_base{
  my $self = shift;
  my $nseq = shift;
  my $pos  = shift;
  if (!defined($pos)){
    $pos = $nseq;
    $nseq = 0;
  }
  my $result = substr($self->{'seqs'}[$nseq]{'seq'}, $pos, 1);
}

sub convert_pos_to_alig{
  my $self = shift;
  my $pos = shift;
  my $name = shift || 0;
  my $result = 0;
  my $counter = 0;
  while ($self->{'seqs'}[$name]{'seq'} =~ /(.)/g && $counter < $pos){
    my $r = $1;
    $counter++ if $r =~ /\w/;
    $result++;
  }
  $result--;
  return $result;
}

sub get_seq{
  my $self = shift;
  my $name = shift;
  my $result = '';
  if (ref($name) eq 'ARRAY'){
    foreach my $el (@$name){
      $result .= $self->get_seq(uc($el));
    }
  }
  else{
    my $done = 0;
    foreach my $prot (keys %{$self->{'byname'}}){
      if ($prot =~ /^$name$/i){
        my $myseq = $self->{'byname'}{$prot};
        my @t = $myseq =~ /(.{0,60})/g;
        $myseq = join("\n", @t);
        $result .= ">$prot\n$myseq\n";
        $done = 1;
      }
    }
    if (!$done){
      if (_isinteger($name)){
        $result = $self->{'seqs'}[$name];
      }
      else{
        warn("Unknown identifier $name\n");
      }
    }
  }
  return $result;
}

sub _isinteger{
  my $n = shift;
  my $result = 1;
  $result = 0 if ($n =~ /\D/);
  return $result;
}

sub prepare_for_alignment{
  my $self = shift;
  my $name = shift;
  my $from = shift;
  my $to = shift;
  my $result = "";
  if (ref($name) eq 'ARRAY'){
    foreach my $el (@$name){
      $result .= $self->prepare_for_alignment(uc($el), $from, $to);
    }
  }
  elsif ($name eq 'all'){
    for (my $i = 0; $i < scalar(@{$self->{'seqs'}}); $i++){
      $result .= $self->prepare_for_alignment($i, $from, $to);
    }
  }
  else{
    my $done = 0;
    foreach my $prot (keys %{$self->{'byname'}}){
      if ($prot =~ /^$name$/i){
        my $myseq = substr($self->{'byname'}{$prot}, $from, $to-$from);
        $result .= "[($prot) ($myseq)]\n";
        $done = 1;
      }
    }
    if (!$done){
      if (_isinteger($name)){
        my $n = $self->{'seqs'}[$name];
        my $myseq = substr($n->{'seq'}, $from, $to-$from);
        $result .= "[($n->{'name'}) ($myseq)]\n";
      }
      else{
        warn("Unknown identifier $name\n");
      }
    }
  }
  return $result;
}

sub seqlen{
  my $self = shift;
  my $nseq = shift || 0;
  my $result = 0;
  if (exists($self->{'seqs'}[$nseq]{'slen'})){
    $result = $self->{'seqs'}[$nseq]{'slen'};
  }
  else{
    $result = length($self->{'seqs'}[$nseq]{'seq'});
    $self->{'seqs'}[$nseq]{'slen'} = $result;
  }
  return $result;
}

sub split_fasta{
  my $self = shift;
  my $ext = shift || '';
  foreach my $seq (@{$self->{'seqs'}}){
    my $name = $seq->{'name'};
    my $s = $seq->{'seq'};
    my ($fname) = $name =~ /(\S+)/;
    $fname .= $ext;
    open (OUT, ">$fname") or die ("Cannot open $fname: $!\n");
    print OUT ">$name\n$s\n";
    close OUT;
  }
}

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

1;


