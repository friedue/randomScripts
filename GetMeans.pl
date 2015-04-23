#! /usr/bin/perl -w

use strict;
use Bio::DB::BigWig 'binMean';
use IO::File;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

my %opt = (
           'bw'    => undef, # big-wig file
           'nb'    => undef, # nb of bins
           'format' => "bed",
	   'summary' => undef,
	   'individual' => undef,
	   'show_zeros' => undef,
	   'verbose'    => undef,
	   'help'    => undef,
          );
GetOptions(
           "bw=s"     => \$opt{bw},
           "nb=i"     => \$opt{nb},
           "format=s" => \$opt{format},
	   "verbose"  => \$opt{verbose},
	   "individual"  => \$opt{individual},
	   "summary"  => \$opt{summary},
	   "show_zeros"  => \$opt{show_zeros},
           "help"     => \$opt{help},
           );


if (!$opt{bw} || $opt{help}) {
    print STDERR "For each chromosomal region (in bed or gff-file) this script calculates the regional mean over scores provided in a bigwig file (-bw)\n";
    print STDERR "Typically this score might be a coverage.\n\n";
    print STDERR "There is an option to split each region into a number of bins (-nb). Then mean score is then calculated for\n";
    print STDERR "each region and every bin.\n";
    print STDERR "If -individual is specified then the script will output the individual Means for each region/bin in the location file. \n";
    print STDERR "If -summary is specified then this script will output the averages over all the regions in the location file. \n";
    print STDERR "If -show_zeros is specified then this script will output 0 for regions which are not covered by the bigwig file over. \n";
    print STDERR "Usage: $0 -bw <big-wig files>  [-nb <number of bins>] [-summary] [-individual] [-format {bed|gff|...}] < locations.bed\n";
    exit(1);
}


my $wig = Bio::DB::BigWig->new(-bigwig=>$opt{bw});


my %average = ();
my %std     = ();
my $N       =  0;


# loop over regions provided as STDIN
while (<>) {
    if ($_ =~ /^\#/) { next; }
    if ($_ =~ /^$/) { next; }
    chomp;
    my %loc = GetLocation($_,\%opt);
    my $length = $loc{end} - $loc{start}  + 1;

    if ($opt{verbose}) { PrintLocation("#inspect region:", \%loc); }
    $N++;

    my $iterator = undef;

    if ($opt{nb}) { 
	###### if regions are to be binned

	my $type="bin:$opt{nb}";
#	$iterator = $wig->features(-iterator=>1, -seq_id=>$loc{chr},-start=>$loc{start},-end=>$loc{end}, -type=>$type);
	$iterator = $wig->features(-iterator=>1, -seq_id=>$loc{chr},-start=>$loc{start},-end=>$loc{end}, -type=>'summary');

	if (!$iterator) { PrintLocation("ERROR. Iterator is empty:", \%loc); exit(1);}

	while (my $it = $iterator->next_seq) {
#		print Dumper($it);
	    
	    my $seqid     = $it->seq_id;
	    my $c_start   = $it->start;
	    my $bin_width = $it->length/$opt{nb};
	    
	    my $stats     = $it->statistical_summary($opt{nb});
	    #print Dumper($stats);

	    if (!$stats) { 
		if (!$opt{individual}) { next ; } 
		
		# the following is only necessary if zero entries are to be output
		# e.g for proper calculation of the median
		my $start = $c_start;
		my $score = 0;
		for (my $jb=0; $jb < $opt{nb}; $jb++) {
		    my $end = $start + $bin_width - 1;
		    print join "\t", ($loc{chr}, int($start), int($end),$score, $jb,"\n"); 
		    $start   += $bin_width;
		} 
		next;
	    } 
	    

	    if ($opt{nb} ne @$stats) {
		PrintLocation("ERROR",\%loc);
		print "ERROR: number of bins (" . @$stats . " is not equal to opt{nb} ($opt{nb})\n";
		exit(1);
	    }

	    my $start = $c_start;
	    my $ib=0;
	    for my $s (@$stats) {
		
		my $jb=$ib;
		if ($loc{strand} eq "-") { $jb   = $opt{nb} - 1 - $ib; }
		
		my $score  = binMean($s);
		if (!$score) { $score = 0; }
		
		$average{$jb}+=$score;
		$std{$jb}+=($score*$score); 
		
		
		if ($opt{individual}) { 
		    my $end   = $start + $bin_width - 1;
		    print join "\t", ($loc{chr}, int($start), int($end),$score, $jb,"\n"); 
		} #bedgraph
		
		$start += $bin_width;
		$ib++;
	    } 
	    
	} # end loop over iterator when binning

	
    } else {

	###### if regions are _not_ to be binned

	$iterator = $wig->features(-iterator=>1, -seq_id=>$loc{chr},-start=>$loc{start},-end=>$loc{end} );
	if (!$iterator) { PrintLocation("ERROR. Iterator is empty:", \%loc); exit(1);}
	#print Dumper($iterator);

	my $found=0;  # keep track of bw-entries found in region %loc
	while (my $it = $iterator->next_seq) {

	    $found++;

	    my $a = $it->start;
	    my $b = $it->end;
	    my $L = $b - $a + 1;
	    my $score=$it->score;
	    if ($opt{individual}) { print join "\t", ($loc{chr}, $a, $b,$score,"\n");  } #bedgraph-output for individual

	      
	    for (my $i=$a; $i<=$b; $i++) {
		
		my $ib   = $i - $loc{start};
		if ($loc{strand} eq "-") { $ib   = $loc{end} - $i; }
		
		$average{$ib}+=$score;
		$std{$ib}    +=($score*$score); 
		#print "$ib $score $average{$ib} $std{$ib}\n";
	    }

	}


	if (!$found) { 
	    # nothing was found in region %loc
	    if (!$opt{individual}) { next ; } 
#	    if ($opt{verbose}) { PrintLocation("#Warning: Region is not found in $opt{bw}", \%loc); }
	    my $score=0;
	    print join "\t", ($loc{chr}, $loc{start}, $loc{end},$score,"\n");  #bedgraph-output for individual
	} 

    } # endif $opt{nb}
    
} # end while (<>)



if ($opt{summary}) {

    if ($N<1) { print STDERR "WARNING! Too few regions to calculate summary. N=$N"; next; }

    foreach my $ib (sort {$a <=> $b} keys %average) {
	my $aver  = 0; if (exists $average{$ib}) { $aver = $average{$ib}/$N; }
	my $error = 0; 
	if (exists $std{$ib}) { 
	    $error = $std{$ib} - $N*$aver*$aver;
	    if ($N>1) { $error /= ($N-1); }
	    if ($error>0) { $error = sqrt($error/$N); }
	}
	
	print join "\t" , ("average ", $ib,$aver,$error,$N,"\n");
    }
}
exit(0);

###########################
sub GetLocation {
    
    my ($line,$opt) = @_;
    my %loc = ();
    
    if ($opt->{format} =~ /gff/i) {
	my @word = split('\t',$line);
	$loc{chr}    = $word[0];
	$loc{start}  = $word[3];
	$loc{end}    = $word[4];
	$loc{strand} = $word[6];
    }

    if ($opt->{format} =~ /bed/i) {
	my @word = split('\t',$line);
	$loc{chr}    = $word[0];
	$loc{start}  = $word[1];
	$loc{end}    = $word[2];
	$loc{strand} = "+";
	if (@word>5) { $loc{strand} = $word[5]; }
    }

    #print join "\t", ($loc{chr},$loc{start}, $loc{strand}), "\n";

    return %loc;
}
    
sub PrintLocation {
    my ($message,$loc) = @_;
    print STDERR join "\t", ($message,$loc->{chr},$loc->{start}, $loc->{end}, $loc->{strand}), "\n";
}
