#!/usr/bin/perl -w
# input: file with chromosome sizes (either .fa or .sizes)
# output: bedGraph file with tiling of chomosomes in intervals of size $length

use strict;
use Getopt::Long;

my $length = undef;
my $shift = undef;
my $help  = undef;
GetOptions("length=i" => \$length,
           "shift=i"  => \$shift,
           "help"     => \$help,
	   );

if ($help || !$length) {
    print STDERR "This script generates a tiling of a genome, whose chromosome sizes are specified in a tab-separated file.\n";
    print STDERR "Usage: $0 -length <int> [-shift <int>] genome.fai\n";
    exit(1);
}

if (!$shift) { $shift=$length; } 

while (<>) {
    my @word = split('\t');
    my $chr = $word[0];
    my $sl  = $word[1];

    print STDERR "$chr $sl\n";
    for (my $start=0; $start <= $sl - $length; $start+=$shift) {
	my $end = $start + $length;
	print join "\t", ($chr,$start,$end,"\n");
    }
}

exit(0);





