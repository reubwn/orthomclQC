#!/usr/bin/perl -w

=head1 TITLE

orthomclFilter_consensusAlignment.pl

=head1 SYNOPSIS

Blasts each n-some OG versus a DB of consensus seqs from all OG alignments
Attempts to correct for the presence of spurious n-some OGs caused by the presence of split-protein alignments

Takes as input:

=over

=item
a) groups.txt file

=item
b) consensusSeqs.fasta

=item
c) the BLASTp DB name (may be the same as above)

=item
d) the OG size class you wish to filter

=back

Returns:

=over

=item
a) correctedGroups.groups -> a new groups file with those OGs inferred to be erroneous REMOVED

=item
b) .descriptions -> blast info for ALL hits

=item
c) .erroneousTwosomes -> blast info for inferred erroneous twosomes

=item
d) .genuineTwosomes -> blast info for those queries that did not have any good hits

=item
e) .log -> more info per hit, plus at the end of this file are lists of both 'erroneous' and 'genuine' n-somes

=back

=head1 USAGE

orthomclFilter_consensusAlignment.pl <groups.txt> <input.fasta> <DB_name> <OG_sizeClass>

=head1 AUTHOR

RW Nowell, November 2013

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Sort::Naturally;
use Bio::Tools::Run::StandAloneBlast;

my %options;
my $results = GetOptions (\%options,
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV == 0) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_consensusAlignment #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

my $groupsFile = $ARGV[0];
my $consensusSeqsFile=$ARGV[1];
my $db = $ARGV[2];
my $OG_sizeClass = $ARGV[3];

## get the consensus of each twosome OG
## (first get OG name of each, extract seq from consensusSeqs)

open(GROUPS, $groupsFile) or die "\n\tCannot open file '$groupsFile': $!\n";
chomp(my @GROUPS=<GROUPS>);
close GROUPS;

print "\tGetting accessions of all $OG_sizeClass-somes...\n";
my %groupsHash;
foreach (@GROUPS) {
	my @a = split(/\: /,$_);
	chomp(my @b = split(/\s+/, $a[1]));
	$groupsHash{$a[0]} = (\@b)  ## k= OG tag, val= array ref of gid|accession's
}

my %n_somesHash;
foreach my $k (keys %groupsHash) {
	my @a = @{$groupsHash{$k}};
	if (@a == $OG_sizeClass) {  ## if OG is of size class x...
		$n_somesHash{$k} = (\@a);  ## ... put into %twosomesHash, where k= OG number, val= array ref of gid|accessions
	}
}
print "\tThere are ".(keys %n_somesHash)." $OG_sizeClass-some OGs found in the file: '".$groupsFile."'\n";

print "\tCreating sequence lookup DB from $ARGV[1]...\n";
my %consensusSeqsHash;
my $in = Bio::SeqIO -> new( -file => $consensusSeqsFile, -format => 'fasta' );

while (my $seq = $in -> next_seq() ) {
  $consensusSeqsHash{( $seq -> display_id() )} = ( $seq -> seq() );  ## k= OG number, val= seq (as string)
}

print "\n\tExtracting ".$OG_sizeClass."-some sequences...\n";
my %x_someSeqsHash;
foreach my $k (keys %n_somesHash) {
	$x_someSeqsHash{$k}=($consensusSeqsHash{$k});  ## k= OG number, v= consensus seq for that twosome
}

print "\tThere are ".(keys %x_someSeqsHash)." $OG_sizeClass-some consensus seqs\n\tWriting to file '$OG_sizeClass-someConsensusSeqs.fasta'...\n";
my $twosomeConsSeqsFile = "$OG_sizeClass-someConsensusSeqs.fasta";

open(SEQS,">".$twosomeConsSeqsFile) or die "\n\tCannot write to file '$twosomeConsSeqsFile': $!\n";
foreach (keys %x_someSeqsHash) {
	my $trimmed = $x_someSeqsHash{$_};
	$trimmed =~ s/^X+|X+$//g;  ## this rather awesome snippet trims leading and trailing 'X' residues from the consensus seqs, to avoid issues with the %query coverage criterion
	print SEQS ">$_\n$trimmed\n";  ## use trimmed seq
}
close SEQS;

print "\n\tPreparing BLASTp analysis... ";
my $outFile=$twosomeConsSeqsFile;
$outFile=~s/.fasta/.BLASTp/;

## invoke tblastn from the command-line **USING THE NEW BLAST+ EXECS**
print "running blast... ";
if (system("blastp", "-query", "$twosomeConsSeqsFile", "-db", "$db", "-evalue", "1e-5", "-outfmt", "5", "-out", "$outFile")!=0) {
	die "\n\tError executing 'tblastn' system command: $!\n";
}
print "finished BLAST.\n\tParsing BLAST output... ";

## now parse blast output
my $search_io = Bio::SearchIO -> new (-format => 'blastxml', -file => $outFile);

my %errorHash;
my %genuineHash;
my $parserOutFile=$twosomeConsSeqsFile;
$parserOutFile=~s/.fasta//;
my $allHits=$parserOutFile.".descriptions";
my $erroneousTwosomes=$parserOutFile.".erroneous_".$OG_sizeClass."-somes";  ## put into this file those queries that DO have a hit that pass the 80/80 rule (ie, putative erroneous twosome OGs)
my $log=$parserOutFile.".log";
open (ALL, ">".$allHits) or die "- Cannot open file: $!\n";
open (ERROR, ">".$erroneousTwosomes) or die "- Cannot open file: $!\n";
open (LOG, ">".$log) or die "- Cannot open file: $!\n";
print ALL "query_name\tq_len\te-value\t%_id\thit_name\thit_len\taln_len\tpercent_query_coverage\n";
print ERROR "query_name\tq_len\te-value\t%_id\thit_name\thit_len\taln_len\tpercent_query_coverage\n";

while( my $result = $search_io->next_result ) {
	while ( my $hit = $result->next_hit ) {
		while ( my $hsp = $hit->next_hsp ) {

			my $percent_id = $hsp->percent_identity;  ## define % identity between query and hit
			my $rounded_id = sprintf("%.2f",$percent_id);

			## defines percent_coverage as the number of residues of the query
			## participating in the hsp / the total query length;
			## ie, over 80% of the query must share over 80% identity with the subject
			## to pass the filter
			my $percent_q_coverage=((($hsp->length('query')/($result->query_length)))*100);
			my $rounded_pqc = sprintf("%.2f",$percent_q_coverage);

			print LOG "- Query = ".$result->query_description."\n";
			print LOG "- Query length = ".$result->query_length."\n";
			print LOG "- Hit = ".$hit->name."\n";
			print LOG "- Length of hit seq = ".$hit->length."\n";
			print LOG "- Hit start -- end position = ".$hsp->start('hit')." -- ".$hsp->end('hit')."\n";

			print ALL $result->query_description,"\t",$result->query_length,"\t",$hsp->evalue,"\t",$rounded_id,"\t",$hit->name,"\t",$hit->length,"\t",$hsp->length('total'),"\t",$rounded_pqc,"\n";

			## if 80/80 rule is passed then hit is inferred to be an ortholog of the query.
			if (($percent_q_coverage>=80) and ($percent_id>=80) and (($result->query_description) ne ($hit->name))) {  ## but discounting hits to self!
				$errorHash{($result->query_description)}=($hit->name); ## k= OG number, val= number of the 'hit' OG
				print ERROR $result->query_description,"\t",$result->query_length,"\t",$hsp->evalue,"\t",$rounded_id,"\t",$hit->name,"\t",$hit->length,"\t",$hsp->length('total'),"\t",$rounded_pqc,"\n";
			}
			print LOG "- Length of HSP = ".$hsp->length('query')."\n";
			print LOG "- %-identity of HSP = ".$hsp->percent_identity."\n";
			print LOG "- Percent query coverage = ".$percent_q_coverage."\n\n";
		}
	}
}
print LOG "\n
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# list of putative error $OG_sizeClass-somes #
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\nquery_OG\thit_OG
";

foreach (nsort keys %errorHash) {
	print LOG "$_\t$errorHash{$_}\n";
}

close ALL;
close ERROR;
close LOG;

print "finished.\n\n\tThere were ".(keys %errorHash)." query seqs that had good hits to another consensus alignment in the dataset\n";
print "\n\tNow removing these OGs from the file ".$groupsFile."... ";

my $correctedGroupsFile = $groupsFile;
$correctedGroupsFile =~ s/.*/correctedFor_$OG_sizeClass-somes.groups/;
open (CORRECTED,">".$correctedGroupsFile) or die "\n\tCannot write to file '$correctedGroupsFile': $!\n";

my %sortableGroupsHash;  ## to sort CORRECTED groups file into numerical order
my $tag;
foreach my $k (nsort keys %groupsHash) {
	# my @a = split(/\_/, $k);
	# $tag = $a[0];  ## puts '27gen' into $tag
	# $sortableGroupsHash{$a[1]} = $groupsHash{$k};  ## key= OG number; val= group

	my @a = @{ $groupsHash{$k} };
	unless (exists $errorHash{$k}) {
		print CORRECTED $k.": @a\n";
	}
}

# print "\tTag is: ".$tag."\n";

# # foreach my $k (sort {$a<=>$b} keys %sortableGroupsHash) {
	# my @a = @{$sortableGroupsHash{$k}};
	# my $search_for = join('_', $tag,$k);  ## because the keys in %errorHash are of type '27gen_1', whereas keys in %sortableGroupsHash are of type '1'
	# unless (exists $errorHash{$search_for}) {  ## if $search_for in hash %errorHash also exists
		# print CORRECTED $tag."_".$k.": @a\n";
	# }
# }
close CORRECTED;

print "done.\n\n\tCorrected groups file is: ".$correctedGroupsFile."\n\n";
