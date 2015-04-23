#! /usr/bin/perl -w


use strict;
use Getopt::Long;

my %opt = (
           'up'        => 0,
	   'down'      => 0,
           'start'     => undef,
	   'end'       => undef,
	   'shift'       => undef,
	   'eps_score'     => undef,
	   'thresh_score'  => undef,
	   'binarize_score'=> undef,
	   'add_score'     => undef,
	   'add_score_ID'  => undef,
	   'transform_score'  => undef,
	   'scale_score'  => undef,
           'c1'     => 0,  # columns in add_score that define region
	   'c2'     => 1,  # columns in add_score that define score
	   'end'       => undef,
           'chr_format'  => undef,
           'help'        => undef,
           'verbose'     => undef,
          );

GetOptions(
           "up=i"   => \$opt{up},
           "down=i" => \$opt{down},
           "start"  => \$opt{start}, # upstream/downstream around start (TSS)
           "end"    => \$opt{end},   # upstream/downstream around end (TES)
           "shift=i"    => \$opt{shift},   # shift both up and down by shift
           "eps_score=f"   => \$opt{eps_score},
           "thresh_score=f"   => \$opt{thresh_score},
           "binarize_score"   => \$opt{binarize_score},
           "add_score=s"      => \$opt{add_score}, # list with scores to added
           "add_score_ID=s"   => \$opt{add_score_ID}, # ID-field to be added
           "transform_score=s"   => \$opt{transform_score}, # transformation for score: log, ...
           "scale_score=f"   => \$opt{scale_score}, # scale factor for score (before all other transformation)
           "c1=i"   => \$opt{c1}, #  columns in add_score that define region
           "c2=i"   => \$opt{c2}, #  columns in add_score that define score
           "help"          => \$opt{help},
           "chr_format=s"  => \$opt{chr_format},
	   "verbose"       => \$opt{verbose}
           );

if ($opt{help}) {
    print STDERR "Usage: $0 -options file.bed\n";
    print STDERR "\tposition modifications: -up <int> -down <int> -shift <int> [-start or -end ]\n";
    print STDERR "\tscore modifications: -eps <double> -thresh <double> -binarize  -add_score <file> [-c1 <int> -c2 <int> columns in add_score file] -add_score_ID <string> -transform_score <string>\n";
    print STDERR "\tchr modifications: -chr_format <ensembl|ucsc> \n";
    exit(1);
}

if ($opt{start} && $opt{end}) {
    print STDERR "ERROR: cannot specify -start and -end simultaneously!\n";
    exit(1);
}


while(<>) {
    if ($_ =~ /^\#/) {print; next;}
    chomp;

    my %feature = ReadFeature($_);

    # change chromosome name (eg. ensembl<-> UCSC)
    $feature{chr} = chr_name($feature{chr},$opt{chr_format});


    # shift
    if ($opt{shift}) { $feature{start} += $opt{shift}; $feature{end} += $opt{shift}; }

    # upstream/downstream
    # define the two positions around which to go upstream/downstream
    # if $opt{start} or $opt{end} are defined than both positions coincide
    # at the start or the end of the current region
    if ($opt{up} != 0 || $opt{down} != 0) {

	my $pos1  = $feature{start};
	my $pos2  = $feature{end};

	if ($feature{strand} eq "-") {
	    $pos1 = $feature{end};
	    $pos2 = $feature{start};
	}

	if ($opt{start}) { $pos2 = $pos1 }
	if ($opt{end})   { $pos1 = $pos2; }

	if ($feature{strand} eq "-") {
	    $feature{start} = $pos2 - $opt{down};
	    $feature{end}   = $pos1 + $opt{up};
	} else {
	    # default: "+" strand
	    $feature{start} = $pos1 - $opt{up};
	    $feature{end}   = $pos2 + $opt{down};
	}
	
	if ($feature{start}<0) {
	    print STDERR "Warning: Skip negative start: id=$feature{id} start=$feature{start} \n";
	    next;
	}
    }
    
    ####################################################################
    # score modifications
    # scale score
    if ($opt{scale_score}) { $feature{score} *= $opt{scale_score}; }

    # bin score
    if ($opt{eps_score}) { $feature{score} = $opt{eps_score} * round($feature{score}/$opt{eps_score});  }

    # threshold score
    if ($opt{thresh_score} && $feature{score} <= $opt{thresh_score} ) { $feature{score} = 0;  }

    # if only region are of interest set all scores above threshold to 1
    # this will have the effect of merging all scores above threshold
    if ($opt{binarize_score} && $feature{score} > $opt{thresh_score} ) { $feature{score} = 1;  }
   

    ####################################################################
    # score transformation
    if ($opt{transform_score} && $opt{transform_score} =~ /^log/i)  { $feature{score} =  log($feature{score}); }
    if ($opt{transform_score} && $opt{transform_score} =~ /^nlog/i) { $feature{score} = -log($feature{score}); }


    PrintFeature(\%feature);

}

exit(0);


sub round {
    my $number = shift;
    return int($number + .5 * ($number <=> 0));
}

sub ReadFeature {
    my $line  = shift;
    my @field = split('\t', $line);

    my %feature = ();

    # bed file
    $feature{chr}   = $field[0];
    $feature{start} = $field[1];
    $feature{end}   = $field[2];
    $feature{id}    = $field[3];
    $feature{score} = $field[4];
    $feature{strand}= $field[5];

    return %feature;
}


sub PrintFeature {
    my $f = shift;
    print join "\t", ($f->{chr}, $f->{start}, $f->{end}, $f->{id}, $f->{score}, $f->{strand},"\n");
}


sub chr_name{
    my $c      = shift;
    my $format = shift;

    if (!$format) { return $c; }

    if ($format =~ /ensembl/i ) {
        $c =~ s/^chr//; $c =~ s/\bM\b/MT/;
    }
    if ($format =~ /ucsc/i ) {
        if ($c !~ /^chr/) { $c =~ s/^MT/M/; $c =~ s/^/chr/; }
    }

    return $c;
}
