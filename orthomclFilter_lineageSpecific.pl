#!/usr/bin/perl -w

=head1 TITLE

orthomclFilter_lineageSpecific.pl

=head1 SYNOPSIS

Uses tBLASTn to query each putative lineage-specific gene versus the DNA contig data for the set of input genomes.
If a good hit is found (>= 80% identity over >= 80% query length), that gene is excluded from the lineage-specifc geneset.
Also excludes if there is no hit found for an input query.

=head1 USAGE

orthomclFilter_lineageSpecific.pl <lineageSpecific.fasta> <BLAST_db_name> <evalue>

=head1 AUTHOR

RW Nowell, October 2013

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use List::Compare;

my %options;
my $results = GetOptions (\%options,
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV == 0) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_lineageSpecific #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n";

print STDERR "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_lineageSpecific #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n";

open (IN, $ARGV[0]) or die "\n\tCannot open file: $ARGV[0]: $!\n\n";

my $numQueries = 0;
while (<IN>) {
	chomp;
	if ($_ =~ m/^>/) { ## fasta header
		$numQueries++;
	}
}
close IN;

my %seqHash; ## hash for rewriting filtered fasta
my $seqIO = Bio::SeqIO -> new( -file => $ARGV[0], -format => 'fasta');
while ( my $seq = $seqIO -> next_seq() ) {
	$seqHash{ $seq -> display_id() } = ( $seq -> seq() );
}

print "\n\tInput file is: ".$ARGV[0]."\n";
print "\n\tNumber of query sequences: ".$numQueries."\n";
print "\tDoing BLAST... ";

my $blastInput = $ARGV[0];
my $blastdb = $ARGV[1];
my $evalue = $ARGV[2];
(my $blastOutput = $ARGV[0]) =~ s/fasta/TBLASTN/;

print "writing BLAST xml to: ".$blastOutput."\n";

## invoke tBLASTn from the command-line
if ( system ("tblastn", "-query", "$blastInput", "-db", "$blastdb", "-evalue", "$evalue", "-outfmt", "5", "-out", "$blastOutput") !=0 ) {
    die "- Error executing system 'tblastn' command: $!\n";
}

print "\n\tFinished BLASTing, now parsing results...\n";

(my $noHits = $ARGV[0]) =~ s/\.fasta/\.noHits/; ## what is this? no hits to any genome, not even self?
(my $oneHits = $ARGV[0]) =~ s/\.fasta/\.oneHit/; ## some details about each hit
(my $fasta = $ARGV[0]) =~ s/\.fasta/\.filtered.fasta/; ## rewrite sequences that pass the filter to another .fasta
(my $log = $ARGV[0]) =~ s/\.fasta/\.log/; ## logfile
my %oneHitHash;
my %hitsHash;
my %noHitHash;

open (FASTA, ">".$fasta) or die "\n\tCannot open file $fasta: $!\n";

my $searchIO = Bio::SearchIO -> new (-format => 'blastxml', -file => $blastOutput);

while( my $result = $searchIO->next_result ) {

	## some queries return no hits - I'm not sure why.. but these are also exluded from final LS lists
	if ($result -> num_hits == 0) {
		$hitsHash{$result->query_description} = ();
		print "\n\tQuery: ".$result->query_description.": No hit found!";
	}

	while ( my $hit = $result->next_hit ) {
		while ( my $hsp = $hit->next_hsp ) {

			## calculate % identity
			my $percent_id = $hsp->percent_identity;
			my $rounded_id = sprintf("%.2f",$percent_id);
			## calculate % query coverage
			my $percent_q_coverage = ((($hsp->length('query')/($result->query_length)))*100);
			my $rounded_coverage = sprintf("%.2f",$percent_q_coverage);
			## get query genome tag
			my $string = ($result->query_description);
			my @a = split(/\|/, $string);
			my $q_gid = $a[0];

			## if a good match is found (>= 80% identity over >=80 query length) in a
			## genome other than the genome from which the sequence originated,
			## the sequence should not be called as lineage-specific:
			if ( ($percent_id >= 80) and ($percent_q_coverage >= 80) and ($q_gid ne $hit->name) ) {
				print "\n\tQuery: ".$result->query_description." has a match in genome: ".$hit->name." with %-identity = ".$rounded_id."%, coverage = ".$rounded_coverage."%, E-value = ".$hsp->evalue;
				$hitsHash{$result->query_description} = $hit->name;
			} else {
				print "\n\tQuery: ".$result->query_description." appears to lineage-specific";
			}
		}
	}
	print "\n\t~~~\n";
}

print "\n\tThese sequences should be demoted from the lineage-specific gene set:\n";

foreach (sort keys %hitsHash) {
	print "\t\t".$_."\n";
}

print "\n\tLeaving these sequences as truly lineage-specific:\n";

## print to filtered .fasta the surviving, truly lineage-specific protein seqs
foreach (sort keys %seqHash) {
	unless (exists $hitsHash{$_}) {
		print "\t\t".$_."\n";
		print FASTA ">".$_."\n".$seqHash{$_}."\n";
	}
}

print "\n\tNumber of inferred LS sequences before filter: ".(keys %seqHash)."\n";
print "\tNumber of inferred LS sequences after filter: ".( (keys %seqHash) - (keys %hitsHash) )."\n\n";

print STDERR "\n\tInput file is: ".$ARGV[0]."\n";
print STDERR "\tNumber of inferred LS sequences before filter: ".(keys %seqHash)."\n";
print STDERR "\tNumber of inferred LS sequences after filter: ".( (keys %seqHash) - (keys %hitsHash) )."\n\n";

close FASTA;

print "\tFinished.\n\n";

__END__
