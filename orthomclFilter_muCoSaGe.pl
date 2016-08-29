#!/usr/bin/perl -w

=head1 TITLE

orthomclFilter_muCoSaGe.pl

=head1 SYNOPSIS

This script actually only counts and prints (to fasta) those groups from the .groups file that
have multiple co-orthologs from the same genome within the same group.

The script orthomclTCoffeeAlign_AA.pl can then be used to create nexus alignments for these groups.

Subsequent visual checking of these alignments can then be done in an alignment viewer such as Geneious.

=head1 USAGE

orthomclFilter_muCoSaGe.pl <groups.txt> <goodProteins.fasta> > muCoSaGe.log

=head1 AUTHOR

RW Nowell, November 2013

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

use Bio::SeqIO;
use Sort::Naturally;

my %options;
my $results = GetOptions (\%options,
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV == 0) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_muCoSaGe #
\t##~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

print STDERR "
\t##~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_muCoSaGe #
\t##~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

open (GROUPS, $ARGV[0]) or die "\n\t$!\n\n";
chomp (my @GROUPS = <GROUPS>);
close GROUPS;

my %GROUPS;
my %groups;
my %goodProteinsHash;

## hash proteins from goodProteins
print STDERR "\tHashing sequences from ".$ARGV[1]."\n";
my $in = Bio::SeqIO -> new ( -file => $ARGV[1], -format => 'fasta' );
while ( my $seq = $in -> next_seq() ) {
	$goodProteinsHash{ $seq->display_id() } = $seq->seq(); ## key= gid|prid, val= seq (as string)
}
print STDERR "\tNumber of proteins: ".(keys %goodProteinsHash)."\n\n";

foreach (@GROUPS) {
	my @a = split (/\: /, $_);
	my @b = split (/\s+/, $a[1]);

	$GROUPS{$a[0]} = \@b; ## key= group ID, val= array of participating proteins

	my @c;
	my %genomes;

	foreach (@b) {
		my @d = split (/\|/, $_);
		push (@c, $d[0]);
	}

	foreach (@c) {
		$genomes{$_}++; ## key= genome ID; val= count of that GID in that group.
	}

	$groups{$a[0]} = \%genomes; ## key= group IG; val= hash as above
}

my %groupsWithParalogs;

foreach my $group (nsort keys %groups) { ## nsort sorts 'naturally', ie lexically then numerically :)
	my %genomes = %{ $groups{$group} };

	foreach (sort { $genomes{$b} <=> $genomes{$a} } keys %genomes) {
		if ($genomes{$_} > 1) { ## ie if there are multiple sequences from the same genome within that group
			$groupsWithParalogs{$group}++; ## key= group ID; val= count of the # genomes that have >1 member in that group
		}
	}
}

## print all groups with muCoSaGe to fasta files
foreach (nsort keys %groupsWithParalogs) {
	my @members = @{ $GROUPS{$_} };

	open (FASTA, ">".$_.".fasta") or die "\n\t$!\n\n";
	foreach (@members) {
		print FASTA"\>".$_."\n".$goodProteinsHash{$_}."\n";
	}
	close FASTA;
}

print "\n\tNumber of groups in total: ".(keys %groups)."\n";
print "\tNumber of groups containing multiple sequences from the same genome: ".(keys %groupsWithParalogs)."\n\n";

print STDERR "\tNumber of groups in total: ".(keys %groups)."\n";
print STDERR "\tNumber of groups containing multiple sequences from the same genome: ".(keys %groupsWithParalogs)."\n\n";

__END__
