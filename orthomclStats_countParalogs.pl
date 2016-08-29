#!/usr/bin/perl -w

=head1 TITLE

orthomclStats_countParalogs.pl

=head1 SYNOPSIS

Simply counts the number of groups that have multiple members from the same genome (ie, apparent paralogs).

=head1 USAGE

orthomclStats_countParalogs.pl <groupsfile>

=head1 AUTHOR

RW Nowell, October 2013

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

use Sort::Naturally;

my %options;
my $results = GetOptions (\%options, 
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV == 0) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclStats_countParalogs #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

print STDERR "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclStats_countParalogs #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

open (GROUPS, $ARGV[0]) or die "\n\t$!\n\n";
chomp (my @GROUPS = <GROUPS>);
close GROUPS;

my %groups;

foreach (@GROUPS) {
	my @a = split (/\: /, $_);
	my @b = split (/\s+/, $a[1]);
	
	my @c;
	my %genomes;
	
	foreach (@b) {
		my @d = split (/\|/, $_);
		push (@c, $d[0]);
	}
	foreach (@c) {
		$genomes{$_}++; ## key= genome ID; val= count of that GID in that group.
	}
	$groups{$a[0]} = \%genomes;
}

my %groupsWithParalogs;

foreach my $group (nsort keys %groups) { ## nsort sorts 'naturally', ie lexically then numerically :)
	my %genomes = %{ $groups{$group} };
	
	my $flag = 0;
	
	foreach (sort { $genomes{$b} <=> $genomes{$a} } keys %genomes) {
		if ($genomes{$_} > 1) {
			print "\tGroup $group has $genomes{$_} sequences for genome $_\n";
			$groupsWithParalogs{$group}++;
			$flag = 1;
		}
	}
	print "\t~~~\n" if $flag == 1;
}

print "\n\tNumber of groups in total: ".(keys %groups)."\n";
print "\tNumber of groups containing multiple sequences from the same genome: ".(keys %groupsWithParalogs)."\n\n";

print STDERR "\tNumber of groups in total: ".(keys %groups)."\n";
print STDERR "\tNumber of groups containing multiple sequences from the same genome: ".(keys %groupsWithParalogs)."\n\n";

__END__