#!/usr/bin/perl -w

=head1 TITLE

orthomclFilter_almostCore.pl

=head1 SYNOPSIS

This script queries those OGs which are 'almost core', in that they are
present in n-1 of the participating genomes.

It does it by taking each 'almost core' group, extracting either the P.s.
phaseolicola 1448A, syringae B728a or tomato DC3000 orthologous protein sequence,
and querying this vs the DNA contig data for the non-participating genome, using tBLASTn.

If a good hit is found in the apparently non-participating genome, that group is
promoted to the fully core group, and is counted towards the core genome of the
input genomes.

A good hit is defined as >80% identity over >80% of the query sequence's length,
or, if the hit regions is very close to the terminus of a contig, >90% sequence identity.

A logfile is written to STDOUT, and a corrected .groups file is also produced,
where any group promoted to core has an additional 'dummy' protein ortholog inserted into it,
called 'name_of_nonparticipating_genome|almostCoreDummySeq' (ie has no sequence information attached).

=head1 DISCLAIMER

Obviously, this script will need to be (fairly extensively) modified depending on the strain names etc.
of the dataset being used.

=head1 USAGE

orthomclFilter_almostCore.pl <.groups file> <goodProteins.fasta> <tblastn evalue> > almostCore.logfile

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
\t##~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_almostCore #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~##\n";

my %groupsHash;
my %genomesHash;
my %goodProteinsHash;
my %hitHash;
my $almostCoreCount = 0;
my $hitCount = 0;
my $close2contigEnd = 0;

open (GROUPS, $ARGV[0]) or die "\n\tCannot open the file '$ARGV[0]': $!\n\n";
while (<GROUPS>) {
	my @a = split(/\: /, $_);
	my @b = split(/\s+/, $a[1]);

	$groupsHash{$a[0]}=(\@b); ## key= group ID; val= array of member proteins

	foreach (@b) {
		chomp;
		my @c = split(/\|/, $_);
		$genomesHash{$c[0]}++; ## key= GID; val= count of GID
	}
}
close GROUPS;

## hash proteins from goodProteins
print "\n\tHashing sequences from ".$ARGV[1]."\n";
my $in = Bio::SeqIO -> new ( -file => $ARGV[1], -format => 'fasta' );
while ( my $seq = $in -> next_seq() ) {
	$goodProteinsHash{ $seq->display_id() } = $seq->seq(); ## key= gid|prid, val= seq (as string)
}
print "\tNumber of proteins: ".(keys %goodProteinsHash)."\n\n";

my $numGenomes = (keys %genomesHash);
my @allGenomes = keys %genomesHash;
print "\tTotal number of genomes: ".$numGenomes."\n\n";

foreach my $group_id (keys %groupsHash) {
	my @b = @{$groupsHash{$group_id}};
	my %nonRedundantGenomeHash; ## core groups need ONE member per genome
	my @currentGenomes;

	foreach (@b) {
		chomp;
		my @c = split(/\|/, $_);
		$nonRedundantGenomeHash{$c[0]}++; ## key= GID; val= count
		push (@currentGenomes, $c[0]);
	}

	## for a group to be almost core, there must be n - 1 UNIQUE identities in @currentGenomes:
	if ( (@currentGenomes == ($numGenomes - 1)) and (keys %nonRedundantGenomeHash == ($numGenomes - 1)) ) {
		print "\tGroup ".$group_id." is an almost core group (has ".($numGenomes-1)." participants)\n";

		## need to identify the missing genome:
		my $lc = List::Compare -> new(\@allGenomes, \@currentGenomes); ## @allGenomes will always be >= @currentGenomes
		my @missingGenome = $lc -> get_unique;
		my $missingGenome = "@missingGenome";
		print "\t\tMissing genome: ".$missingGenome."\n";
		my $query;

		## if missingGenome is in PG3 (or if syriB728a is the non-participating genome), use phas1448A ortholog:
		if (($missingGenome =~ m/rhap4220/) or ($missingGenome =~ m/myri2897/) or ($missingGenome =~ m/mors2341/) or ($missingGenome =~ m/mors5269/) or
		    ($missingGenome =~ m/erio2343/) or ($missingGenome =~ m/dend3226/) or ($missingGenome =~ m/neri5067/) or ($missingGenome =~ m/sava3335/) or
		    ($missingGenome =~ m/frax5062/) or ($missingGenome =~ m/daph4219/) or ($missingGenome =~ m/cast4217/) or ($missingGenome =~ m/aesc2329/) or
		    ($missingGenome =~ m/aesc2336/) or ($missingGenome =~ m/aesc2306/) or ($missingGenome =~ m/aesc2250/) or ($missingGenome =~ m/aesc3681/) or
		    ($missingGenome =~ m/aesc2279/) or ($missingGenome =~ m/aesc2113/) or ($missingGenome =~ m/aesc2315/) or ($missingGenome =~ m/ulmi1407/) or
		    ($missingGenome =~ m/cera6109/) or ($missingGenome =~ m/glycR4/) or ($missingGenome =~ m/glycB076/) or ($missingGenome =~ m/taba11528/) or
		    ($missingGenome =~ m/lach301315/) or ($missingGenome =~ m/brou5140/) or ($missingGenome =~ m/mori301020/) or ($missingGenome =~ m/syriB728a/) or
		    ($missingGenome =~ m/tomaDC3000/)) {
			foreach (@b) {
				if ($_ =~ m/phas1448A/) {
					$query = $_;
					print "\t\tMissing genome is from PG3 - Using phaseolicola 1448A ortholog: \n\t\t\t".$query."\n";

					## print to a temp fasta for input into BLAST
					open (TEMPSEQ, ">query.temp.fasta") or die "\n\tCannot write to file 'query.fasta': $!\n\n";
					print TEMPSEQ ">$_\n$goodProteinsHash{$_}";
					close TEMPSEQ;
				}
			}
		## if missingGenome is from PG2 (or if phas1448A is the non-participating genome), use syriB728a ortholog:
		} elsif (($missingGenome =~ m/BRIP34881/) or ($missingGenome =~ m/BRIP34876/) or ($missingGenome =~ m/japo301072/) or ($missingGenome =~ m/pani2367/) or
		    ($missingGenome =~ m/syriFF5/) or ($missingGenome =~ m/apta50252/) or ($missingGenome =~ m/pisi1704B/) or ($missingGenome =~ m/avel013/) or
		    ($missingGenome =~ m/avel037/) or ($missingGenome =~ m/syri7872/) or ($missingGenome =~ m/acer302273/) or ($missingGenome =~ m/syri2339/) or
		    ($missingGenome =~ m/syri2340/) or ($missingGenome =~ m/syri7924/) or ($missingGenome =~ m/cit7/) or ($missingGenome =~ m/BRIP39023/) or
		    ($missingGenome =~ m/papu1754/) or ($missingGenome =~ m/syri642/) or ($missingGenome =~ m/phas1448A/)) {
			foreach (@b) {
				if ($_ =~ m/syriB728a/) {
					$query = $_;
					print "\t\tMissing genome is from PG2 - Using syringae B728a ortholog: \n\t\t\t".$query."\n";
					open (TEMPSEQ, ">query.temp.fasta") or die "\n\tCannot write to file 'query.fasta': $!\n\n";
					print TEMPSEQ ">$_\n$goodProteinsHash{$_}";
					close TEMPSEQ;
				}
			}
		} elsif (($missingGenome =~ m/tomaT1/) or ($missingGenome =~ m/tomaK40/) or ($missingGenome =~ m/tomaMax13/) or ($missingGenome =~ m/toma1108/) or
		    ($missingGenome =~ m/lach302287/) or ($missingGenome =~ m/avii3846/) or ($missingGenome =~ m/actn3739/) or ($missingGenome =~ m/actn3871/) or
		    ($missingGenome =~ m/actnCRAFRU843/) or ($missingGenome =~ m/mors302280/) or ($missingGenome =~ m/mors5261/) or ($missingGenome =~ m/avel631/) or
		    ($missingGenome =~ m/actn302091/) or ($missingGenome =~ m/oryz1_6/) or ($missingGenome =~ m/macuES4326/)) {
			foreach (@b) {
				if ($_ =~ m/tomaDC3000/) {
					$query = $_;
					print "\t\tMissing genome is from PG1 - Using tomato DC3000 ortholog: \n\t\t\t".$query."\n";
					open (TEMPSEQ, ">query.temp.fasta") or die "\n\tCannot write to file 'query.fasta': $!\n\n";
					print TEMPSEQ ">$_\n$goodProteinsHash{$_}";
					close TEMPSEQ;
				}
			}
		}

		## run tBLASTn on the phaseolicola sequence versus the appropriate genome:
		my $evalue = $ARGV[2];
		# my $outname = $group_id.".TBLASTN";

		print "\t\tRunning tBLASTn with following parameters:\n";
		print "\t\t\tQuery name: ".$query."\n";
		print "\t\t\tDB name: ".$missingGenome."\n";
		print "\t\t\tEvalue: ".$ARGV[2]."\n";
		# print "\t\t\tOutfile name: ".$outname."\n";

		if (system("tblastn", "-query", "query.temp.fasta", "-db", "$missingGenome", "-evalue", "$evalue", "-max_target_seqs", "1", "-outfmt", "5", "-out", "query.temp.TBLASTN")!=0) {
			die "\tError executing system command: tblastn\n\n";
		}
		print "\t\tDone BLASTing, now parsing results:\n";

		## parse BLAST output
		my $searchIO = Bio::SearchIO -> new (-format => 'blastxml', -file => 'query.temp.TBLASTN');

		while( my $result = $searchIO->next_result ) {
				while ( my $hit = $result->next_hit ) {
					while ( my $hsp = $hit->next_hsp ) {

					## calculate % identity
					my $percent_id = $hsp->percent_identity;
					my $rounded_id = sprintf("%.2f",$percent_id);
					## calculate % query coverage
					my $percent_q_coverage = ((($hsp->length('query')/($result->query_length)))*100);
					my $rounded_q_coverage = sprintf("%.2f",$percent_q_coverage);

					print "\t\t\tHit name: ".$hit->name."\n";
					print "\t\t\t%-identity: ".$rounded_id."\n";
					print "\t\t\t%-query coverage: ".$rounded_q_coverage."\n";

					## if 80/80 rule is passed then hit is inferred to be an ortholog of the
					## query. HOWEVER, the 80/80 rule can be violated if the hit maps to a
					## truncated protein at the end of a contig - thus, if the hit start or end
					## is within 150bp of a contig end, but the hit itself still has excellent (>90)
					## homology to query, it is assumed to be an ortholog and is counted as a member
					## of the core genome. This kind of thing may be the #1 reason for not being
					## annotated as a gene in the first place.
					if ( ($rounded_id >= 80) and ($percent_q_coverage >= 80) ) {
						print "\t\t* Hit to non-participating genome found!\n\t\tPasses 80/80 rule, OG $group_id inferred to be part of the core genome.\n";
						$hitHash{$group_id} = $missingGenome; ## key= group_id; val= identity of the missingGenome
						$hitCount++;
						last;
					} elsif ( ($hsp->start('hit')<=150) or (($hsp->start('hit')>=(($hit->length)-150))) or ($hsp->end('hit')<=150) or (($hsp->end('hit')>=(($hit->length)-150))) and ($percent_id>=90) ) { ## ie, if the hit start/end is within 150 bp of the beginning/end of the hit contig
						print "\t\t** Hit is found to be close to contig end! Percent-identity of hit >= 90%, OG $group_id inferred to be part of the core genome.\n";
						$hitHash{$group_id} = $missingGenome; ## key= group_id; val= identity of the missingGenome
						$hitCount++;
						$close2contigEnd++;
						last;
					} else {
						print "\t\tNo hit found, OG appears to be genuinely missing from non-participating genome.\n";
						last;
					}
				}
			}
		}
		$almostCoreCount++;
		print "\n\t~~~\n\n";
	} elsif ( (@currentGenomes == $numGenomes) and (keys %nonRedundantGenomeHash == $numGenomes) ) {
			## TODO: calculate the number of original core groups...
	}
}

print "\n\t~~~\n\tNumber of almost core groups: ".$almostCoreCount."\n";
print "\tTotal number of almost core groups that should be promoted to 'core': ".$hitCount."\n";
print "\t[Number of almost core groups with <80% coverage, but >90% identity AND within 150bp of a contig terminus: ".$close2contigEnd."]\n\t~~~\n\n\tThe following groups need adjusting:\n";
foreach (sort {$hitHash{$a} cmp $hitHash{$b}} keys %hitHash) {
	print "\t\t".$_."\t".$hitHash{$_}."\n";
}
print "\n\n\tFinished";

## TODO: rewrite the .groups file based on these corrections - anything with a good hit needs a 'dummy seq' integrated into it's group profile...

## sort %$groupsHash to reflect original layout
my $groupsTag;
my %numberedGroupsHash;
foreach my $group_id (keys %groupsHash) {
	my @a = split(/_/, $group_id);
	$groupsTag = $a[0]; ## 63gen
	$numberedGroupsHash{$a[1]}=$groupsHash{$group_id}; ## key= number; val= array ref in original groupsHash
}

## rewrite the .groups file based on these corrections - anything with a good hit needs a 'dummy seq' integrated into it's group profile...
(my $newGroupsFile = $ARGV[0]) =~ s/.groups/.almostCore.groups/;
open (NEWGROUPSFILE, ">".$newGroupsFile);

foreach my $groupNumber (sort {$a<=>$b} keys %numberedGroupsHash) {
	my $group_id = $groupsTag."_".$groupNumber;
	my @members = sort {$a cmp $b} @{ $numberedGroupsHash{$groupNumber} };

	## if this group is an almost core group needing a correction:
	if (exists $hitHash{$group_id}) {
		print NEWGROUPSFILE $group_id.": @members"." ".$hitHash{$group_id}."|almostCoreDummySeq\n";
	} else {
		print NEWGROUPSFILE $group_id.": @members\n";
	}
}

close NEWGROUPSFILE;

__END__
