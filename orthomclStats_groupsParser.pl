#!/usr/bin/perl -w

=head1 TITLE

orthomclStats_groupsParser.pl

=head1 SYNOPSIS

Decomposes the output from orthMCL (eg 63gen.groups) into lists of
core (1:1:1 orthology) and strain-specific proteins.

Takes the .groups file as the first argument and goodProteins.fasta as the second.

=head1 USAGE

orthomclStats_groupsParser.pl <orthomcl.groups> <goodProteins.fasta>

=head1 AUTHOR

RW Nowell, July 2013

=cut

use strict;

use POSIX;
use Getopt::Long;
use Pod::Usage;
use Sort::Naturally;

my %options;
my $results = GetOptions (\%options,
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV != 2) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclStats_groupsParser #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

my $orthomclInputFile=$ARGV[0];
my $goodProteinsFile=$ARGV[1];
my $outFile=$orthomclInputFile;
$outFile =~ s/.txt//;

open (GP, $goodProteinsFile) or die "\n\t- Cannot open '$goodProteinsFile': $!\n";
my @gp_in = <GP>;
close GP;

my @wholeProteomicDataset;
my %gid_hash;
foreach (@gp_in) {
	chomp;
	if ($_ =~ /^>/) {
		my $header = $_;
		$header =~ s/>//;
		push (@wholeProteomicDataset, $header);

		my @a=split(/\|/,$_);
		my $gid=$a[0];
		$gid=~s/>//;
		$gid_hash{$gid}=();
	}
}

open (MCL, $orthomclInputFile) or die "\n\t- Cannot open '$orthomclInputFile': $!\n";
my @mcl_in = <MCL>;
close MCL;

my @wholeOrthoMCLDataset;
my @OG_size;
foreach (@mcl_in) {
	chomp;
	my @OG = split(/\s+/, $_);
	shift @OG;
	push(@wholeOrthoMCLDataset, @OG);

	my $size = scalar(@OG);
	push (@OG_size, $size);
}

print "
\tTotal number of input genomes: ".(keys %gid_hash)."\n";
foreach (sort {$a cmp $b} keys %gid_hash) { print "\t\t$_\n" };

print "
\tTotal number of proteins across all input genomes: ".@wholeProteomicDataset."
\tNumber of proteins that form an orthologous group: ".@wholeOrthoMCLDataset."
\tNumber of orthologous groups: ".@mcl_in."
";

## get strain specific proteins ##

my %seen;
@seen{@wholeOrthoMCLDataset} = ();

my @strainSpecific;
foreach (@wholeProteomicDataset) {
  push ( @strainSpecific, $_ ) unless exists $seen{$_};
}

print "\n\tTotal number of strain specific proteins: ".@strainSpecific."\n";

my $file_out = "$outFile.strainSpecific";
if (-e $file_out) {
    print "\n\tThe file '$file_out' already exists, would you like to overwrite it? [y/n]: "; ## asks to overwrite
    chomp (my $answer = <STDIN>);
    if (($answer eq "y") or ($answer eq "Y")) {
        unlink $file_out;
    } else {
        die "- Exiting script\n";
    }
}
open (SS, ">$file_out") or die "- $!\n";
foreach (nsort @strainSpecific) {
	print SS "$_\n";
}

## get Core proteins ##

my %result;
my $counter;

foreach (@mcl_in) {
	my @a=split(/\s+/,$_);
	my $group=shift @a;  ## gets OG number tag
	$group=~s/://;
	my @gid_array;
	my @prid_array;

	foreach (@a) {
		my @b = split (/\|/,$_);
		chomp (my $gid = $b[0]);
		chomp (my $prid = $b[1]);
		push (@gid_array, $gid);  ## set up array of genome ids...
		push (@prid_array, $prid);  ## and protein equivalent
	}

	my %query;
	@query{@gid_array}=();  ## sets up query hash

	## asks 'if keys in %wanted equal keys in %query, AND there are NO (==0) keys in %wanted that exist in %query', do...
	## ie checks if keys in %wanted and %query are the same
	if ( @a == (keys %gid_hash) ) {
		if ( (keys %gid_hash == keys %query) and (0 == grep {!exists $query{$_}} keys %gid_hash) ) {
			my $r="@a";
			$result{$group}=$r;  ## pushes to %result
			$counter++;
		}
	}
}

print "\n\tNumber of proteins inferred to be core: ".$counter."\n";

$file_out = "$outFile.core";
if (-e $file_out) {
    print "\n\tThe file '$file_out' already exists, would you like to overwrite it? [y/n]: ";
    chomp (my $answer = <STDIN>);
    if (($answer eq "y") or ($answer eq "Y")) {
        unlink $file_out;
    } else {
        die "- Exiting script\n";
    }
}
open (CORE, ">$file_out") or die "- $!\n";
foreach my $k (nsort keys %result) {
	print CORE "$k: $result{$k}\n";
}

## stats ##

my %size;
foreach (@OG_size) {
	$size{$_}++;
}

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
my $time_now = localtime;

$file_out = "$outFile.stats";
if (-e $file_out) {
    print "\tThe file '$file_out' already exists, would you like to overwrite it? [y/n]: ";
    chomp (my $answer = <STDIN>);
    if (($answer eq "y") or ($answer eq "Y")) {
        unlink $file_out;
    } else {
        die "- Exiting script\n";
    }
}
open (STATS, ">$file_out") or die "- $!\n";

print STATS "
------------------------------------
FREQUENCY DISTRIBUTION of size of OG
OG_size = number of proteins within a given OG
count = frequency of occurence of a given OG size
prop_freq = proportion relative to total number of OG sizes
-----------------------------------------------------------

OG_size\tcount\tprop_freq
";

foreach ( sort {$a <=> $b} keys %size ) {
	my $prop_freq = (($size{$_})/(@mcl_in));
	my $rounded_prop_freq = sprintf("%.4f", $prop_freq);
	print STATS "$_\t$size{$_}\t$rounded_prop_freq\n";
}

my @SSgids;
foreach (@strainSpecific) {
	chomp;
	my @a=split(/\|/, $_);
	my $gid=$a[0];
	push(@SSgids, $gid);
}

my %count;
foreach (@SSgids) {
	$count{$_}++;
}

print STATS "
---------------------------------
NUMBER STRAIN-SPECIFIC PER GENOME
Counts of numbers of inferred strain-specific proteins per input genome
(Before filtering)
-----------------------------------------------------------------------

genome\tcount
";
foreach my $k (sort keys %count) {
	print STATS "$k\t$count{$k}\n";
}

print STATS "

- Input file: $orthomclInputFile
- File created on: $time_now

- Total number of proteins across all input genomes: ".@wholeProteomicDataset."
- Number of proteins that form an orthologous group: ".@wholeOrthoMCLDataset."
- Total number of strain specific proteins: ".@strainSpecific."
- Number of proteins inferred to be core: ".$counter."
";

print "\n\tFinished\n\n";

__END__
