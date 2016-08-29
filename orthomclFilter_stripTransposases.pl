#!/usr/bin/perl -w

=head1 TITLE

orthomclFilter_stripTransposases.pl

=head1 SYNOPSIS

Uses BLAST to identify which proteins are IS elements, transposases or related
mobilome sequences, and removes the OG in which these sequences are found from the
.groups file.

=head1 USAGE

orthomclFilter_stripTransposases.pl <groups.txt> <goodProteins.fasta> <blast DB> <blast evalue> > stripTransposases.logfile

=head1 AUTHOR

RW Nowell, October 2013

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my %options;
my $results = GetOptions (\%options,
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV == 0) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_stripTransposases #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n";

my $groupsFile = $ARGV[0];
my $goodProteinsFile = $ARGV[1];
my $blastDB = $ARGV[2];
my $evalue = $ARGV[3];
(my $blastOutfile = $goodProteinsFile) =~ s/.fasta/.BLASTP/;
(my $newGroupsFile = $groupsFile) =~ s/.groups/.stripTransposases.groups/;
my %hitsHash;

## run tBLASTn on the phaseolicola sequence versus the appropriate genome:
print "\n\tRunning BLASTp with following parameters:\n";
print "\t\tDB name: ".$blastDB."\n";
print "\t\tEvalue: ".$evalue."\n";

if (system("blastp", "-query", "$goodProteinsFile", "-db", "$blastDB", "-evalue", "$evalue", "-max_target_seqs", "1", "-outfmt", "5", "-out", "$blastOutfile")!=0) {
	die "\tError executing system command: blastp\n\n";
}
print "\tDone BLASTing, now parsing results:\n";

## parse BLAST output
my $searchIO = Bio::SearchIO -> new (-format => 'blastxml', -file => $blastOutfile);

while( my $result = $searchIO->next_result ) {
	while ( my $hit = $result->next_hit ) {
		while ( my $hsp = $hit->next_hsp ) {

			## calculate % identity
			my $percent_id = $hsp->percent_identity;
			my $rounded_id = sprintf("%.2f",$percent_id);
			## calculate % query coverage
			my $percent_q_coverage = ((($hsp->length('query')/($result->query_length)))*100);
			my $rounded_q_coverage = sprintf("%.2f",$percent_q_coverage);

			## if 80/80 rule is passed then hit is inferred to be an ortholog of the query.
			if ( ($percent_id >= 80) and ($percent_q_coverage >= 80) ) {
				print "\t~~~\n";
				print "\t\tQuery name: ".$result->query_description."\n";
				print "\t\tHit name: ".$hit->name."\n";
				print "\t\t%-identity: ".$rounded_id."\n";
				print "\t\t%-query coverage: ".$rounded_q_coverage."\n";

				$hitsHash{$result->query_description} = ($hit->name); ## key= accession number of protien inferred to be an IS element etc; val= name of hit.

				last;
			}
		}
	}
}
print "\n\tNumber of IS-related proteins in ".$goodProteinsFile.": ".(keys %hitsHash)."\n\n";

## now, need to find all the group_ids that contain any sequence in %hitsHash
open (GROUPS, $groupsFile) or die "\n\tCannot open $groupsFile: $!\n\n";
my @GROUPS = <GROUPS>;
close GROUPS;

open (NEWGR, ">".$newGroupsFile) or die "\n\tCannot write to $newGroupsFile: $!\n\n";

my %excludeGroups;
my $count = 0;

foreach (@GROUPS) {
	my @a = split (/\: /, $_);
	my @members = split (/\s+/,$a[1]);

	foreach (@members) {
		if ( exists $hitsHash{$_} ) {
			print "\tGroup ".$a[0]." contains the sequence ".$_.", inferred to be ".$hitsHash{$_}."\n";
			$excludeGroups{$a[0]}++;
		}
	}
}
print "\n\tNumber of groups to be excluded: ".(keys %excludeGroups)."\n\n";

foreach (@GROUPS) {
	my @a = split (/\: /, $_);
	unless (exists $excludeGroups{$a[0]}) {
		print NEWGR $_;
		$count++;
	}
}
print "\tNumber of groups in original file '$groupsFile': ".@GROUPS."\n";
print "\tNumber of groups written to '$newGroupsFile': $count\n\n";
print "\n\tFinished.\n\n";

__END__
