#! /usr/bin/perl -w


# This scripts cleans bed-files 
use strict;
use Getopt::Long;

my %opt = (
           'size' => undef,
           'unique_coord' => undef,
           'help'        => undef,
          );

GetOptions(
           "size"   => \$opt{size},
           "unique_coord"   => \$opt{unique_coord},
           "help"          => \$opt{help},
           );

if ($opt{help}) {
    print STDERR "Usage: $0 -options file.bed\n";
    print STDERR "\t-unique_coord : remove redundant coordinates from file\n";
    exit(1);
}


my %sizes=();
if ($opt{size} && -e $opt{size}) {
    open(SIZES, "<$opt{sizes}") || die "Cannot open file $opt{sizes} \n";
    while(<SIZES>) {
	if ($_ =~ /^\#/) {print; next;}
	chomp;
	my @word=split('\s');
	$sizes{$word[0]}=$word[1];
    }

}

my %loc=();
while(<>) {
    if ($_ =~ /^\#/) {print; next;}
    chomp;
    my %feature = ReadFeature($_);

    my $l=$feature{chr} . ":" . $feature{start} . "-" . $feature{end};

    if (exists $loc{$l} && $opt{unique_coord}) { print STDERR "not_unique: $l\n"; next; }

    if (exists $sizes{$feature{chr}}  && $feature{end}>$sizes{$feature{chr}} ) { print STDERR "skip too long: $l\n"; next; }
    
    PrintFeature(\%feature);
    $loc{$l}++;
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
