#!/usr/bin/perl -w

=head1 TITLE

orthomclFilter_createConsensusDB.pl

=head1 SYNOPSIS

Takes the groups file and a) gets the protein seq data for each OG, b) makes an alignment,
c) gets the consensus seq from the alignment, d) saves a flat-file of these consensus seqs, to be subsequently made into a BLAST db

Takes as input: a) the groups.txt file, b) the goodProteins.fasta file, c) an integer for the desired percent
threshold for the consensus alignment (eg 50% consensus threshold, or set to 0% for majority consensus (recommended)).

=head1 WARNINGS

TCoffee might throw a non-fatal warning like:

WARNING: Sequence Seq63 does not contain any residue: automatically removed from the set [WARNING:T-COFFEE]

This is due to the inclusion of 'dummy' seqs after the almostCore filtering step - you can ignore these warnings!

=head1 USAGE

orthomclFilter_createConsensusDB <groups.txt> <goodProteins.fasta> <consensus_threshold> > createConsensusDB.log

=head1 AUTHOR

RW Nowell, November 2013

=cut

use strict;
use Getopt::Long;
use Pod::Usage;

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Sort::Naturally;
use Bio::Tools::Run::Alignment::TCoffee;

my %options;
my $results = GetOptions (\%options,
    'help|h'
    ) or pod2usage(-verbose => 3);

pod2usage(-verbose=>3) if ( ($options{'help'}) or (@ARGV == 0) );

print "
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
\t# orthomclFilter_createConsensusDB #
\t##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##\n\n";

my $groupsFile=$ARGV[0];
my $goodProteinsFile=$ARGV[1];

open(GROUPS, $groupsFile) or die "\n\tCannot open file '$groupsFile': $!\n";
chomp(my @GROUPS=<GROUPS>);
close GROUPS;

print "\tGetting accessions...\n";
my %groupsHash;
foreach (@GROUPS) {
	chomp;
	my @a = split(/\: /, $_);  ## the space after the ':' is important!
	my @b = split(/\s+/, $a[1]); ## splits the gid|accession per line, puts into @b

	$groupsHash{$a[0]} = (\@b);  ## k= OG number (eg 27gen_100 etc), v= ref to array of accessions
}

print "\tCreating sequence lookup DB from $ARGV[1]...\n";
my %goodProteins_hash;
my $in = Bio::SeqIO -> new( -file   => $goodProteinsFile,
                            -format => 'Fasta');

while (my $seq = $in -> next_seq() ) {
  $goodProteins_hash{( $seq->display_id() )} = ( $seq->seq() ); ## key= seq name (gid|prid); val= seq (as string)
}

print "\tGetting seq data per OG...\n";
my %seqsHash;
foreach my $k (nsort keys %groupsHash) {
	my @a = @{$groupsHash{$k}};  ## @a is an array of accessions
	my @seqList;
	foreach (@a) {
		push(@seqList, $goodProteins_hash{$_});
	}
	$seqsHash{$k} = (\@seqList);  ## k= OG number, v= array ref to seqs
}

## need to convert \@seqList from an array of strings -> array of Bio::Seq objects
my %seqsHash2;
foreach my $k (nsort keys %seqsHash) {
	my @seqList_asString=@{$seqsHash{$k}};
	my @seqList_asObjs;

	foreach (@seqList_asString) {
		my $seqobj = Bio::Seq->new( -seq => $_);
		push (@seqList_asObjs, $seqobj);
	}
	$seqsHash2{$k} = (\@seqList_asObjs);  ## k= OG number, v= array ref to array of Bio::Seq objs
}

my $total=(keys %seqsHash2);
my $counter=1;

my @params = ('ktuple' => 1, 'matrix' => 'pam', 'quiet' => '/null/dev');  ## TCoffee params
my $factory = Bio::Tools::Run::Alignment::TCoffee -> new(@params);  ## open new TCoffee factory object

my %alnsHash;
foreach my $k (nsort keys %seqsHash2) {  ## TCoffee takes an array ref to an array of seqs... value of %seqs IS exactly this...
	my $aln = $factory->align($seqsHash2{$k});  ## aligns seqs, happy days

	my $complete = (($counter/$total)*100);
	my $rounded=sprintf("%.2f", $complete);
	print "\tAligning OG '$k', percent complete = $rounded\n";  ## some completedness stats
	$counter++;

	$alnsHash{$k}=$aln;  ## k= OG number, v= Bio::Align::AlignI object
}

print "\n\tDone aligning\n\tThere are ".(keys %alnsHash)." alignments in total\n\tGetting consensus seqs from alignments...\n";

## now get consensus from aln objs in %alnsHash
my $threshold = $ARGV[2];
print "\n\tConsensus threshold set to: $threshold percent\n";
my %consensusHash;
foreach my $k (keys %alnsHash) {
	my $aln = $alnsHash{$k};
	my $consensus = $aln->consensus_string($threshold);  ## get consensus string
	$consensusHash{$k} = $consensus;  ## k= OG name, val= consensus of alignment as string
}

print "\n\tThere are ".(keys %consensusHash)." consensus seqs\n\tWriting consensus seqs to file...\n";

(my $consFile = $groupsFile) =~ s/.groups/.consensusSeqsDB.fasta/;

open (CONSENSUS, ">".$consFile) or die "\n\tCannot write to the file '$consFile': $!\n";
foreach my $k (keys %consensusHash) {
	print CONSENSUS ">".$k."\n".$consensusHash{$k}."\n";
}
close CONSENSUS;

print "\n\tFinished.\n\n";

__END__
