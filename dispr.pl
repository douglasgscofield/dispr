#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;
use Bio::Tools::SeqPattern;

#"class1:F:CTNCAYVARCCYATGTAYYWYTTBYT"
#"class1:R:GTYYTNACHCYRTAVAYRATRGGRTT"

my @o_pf = qw/CTNCAYVARCCYATGTAYYWYTTBYT/;
my @o_pr = qw/GTYYTNACHCYRTAVAYRATRGGRTT/;
my $o_orientation = "FR";
my $o_multiplex = 1;
my $o_min = 300;
my $o_max = 1000;
my $o_bed = "out.bed";
my $o_seq = "out.fa";
my $o_ref = "test.fa";


my $usage = "
$0  --forward FORWARD_PRIMER  --reverse REVERSE_PRIMER --orientation FR|RF|FF|ALL  --reference FASTA_FILE

Applies forward and reverse primer pairs to the reference file,
identifying amplicons within given dimensions.  The primer sequences
are sought paired in specified orientations within a range given
by --min and --max.  The amplicon is measured from the outer extent
of each primer sequence, so includes the primer sequence length.

    --pf FORWARD_PRIMER Forward primer sequence (may be specified more than once)
    --pr REVERSE_PRIMER Reverse primer sequence (may be specified more than once)
    --primers FILE        File containing primer pair(s), as forward/reverse
    --multiplex           If more than one primer pair presented, consider amplicons
                          produced by any possible primer pair (DEFAULT)
    --no-multiplex        If more than one primer pair presented, only consider
                          amplicons produced by each primer pair
    --orientation FR      Orientation of primers, only FR supported for now
    --min bp              Minimum accepted size of amplicon
    --max bp              Maximum accepted size of amplicon

    --both                Orientation of the reference sequence
    --forward
    --reverse

    --ref INPUT_FASTA     Input, Fasta reference sequence in which to find amplicons

    --bed OUTPUT_BED      Output, BED file containing identified amplicon positions
    --seq OUTPUT_FASTA    Output, Fasta sequences containing identified amplicon sequences
    --verbose             Describe actions

";

die "only one primer pair currently supported" if @o_pf > 1 or @o_pr > 1;
die "only FR orientation currently supported" if $o_orientation ne "FR";

sub expand_dot($) {
    # would rather be explicit about '.'
    my $pat = shift;
    $pat =~ s/\./[ACTG]/g;
    return $pat;
}

# http://stackoverflow.com/questions/87380/how-can-i-find-the-location-of-a-regex-match-in-perl
sub match_positions($$) {
    my ($pat, $seq) = @_;
    my @ans;
    while ($$seq =~ /$pat/ig) {
        my ($beg, $end) = ($-[0], $+[0]);
        my $hit = substr($$seq, $beg, $end - $beg);
        print STDERR "match_positions: beg: $beg   end: $end   hit: $hit\n";
        push @ans, [$beg, $end, $hit];
    }
    return @ans;
}

my $seq_pf = Bio::Seq->new(-seq => $o_pf[0], -alphabet => 'dna');
my $seq_pr = Bio::Seq->new(-seq => $o_pr[0], -alphabet => 'dna');
my $sp_pf = Bio::Tools::SeqPattern->new(-seq => $seq_pf->seq(), -type => 'dna');
my $sp_pr = Bio::Tools::SeqPattern->new(-seq => $seq_pr->seq(), -type => 'dna');
my $pat_pff = expand_dot($sp_pf->expand());
my $pat_pfr = expand_dot($sp_pf->revcom(1)->expand());
my $pat_prf = expand_dot($sp_pr->expand());
my $pat_prr = expand_dot($sp_pr->revcom(1)->expand());
$pat_pff = qr/$pat_pff/i;
$pat_pfr = qr/$pat_pfr/i;
$pat_prf = qr/$pat_prf/i;
$pat_prr = qr/$pat_prr/i;
print STDERR "pat_pff: $pat_pff   pat_pfr: $pat_pfr\n";
print STDERR "pat_prf: $pat_prf   pat_prr: $pat_prr\n";

my $stream_pf = Bio::Tools::IUPAC->new(-seq => $seq_pf);
my $stream_pr = Bio::Tools::IUPAC->new(-seq => $seq_pr);
my $fnumber = 0;
++$fnumber while $stream_pf->next_seq();
print STDERR "f: $fnumber expanded sequences\n";
my $rnumber = 0;
++$rnumber while $stream_pr->next_seq();
print STDERR "r: $rnumber expanded sequences\n";


my $in = Bio::SeqIO->new(-file => "<$o_ref", -format => 'fasta');
my $inseq = $in->next_seq();

print STDERR "looking for '$pat_pff' in '".$inseq->display_id()."'\n";

my @hits = match_positions($pat_pff, \$inseq->seq);
#my @hits = match_positions("AACCCTACCTAAACCTCA", \$inseq->seq);
#my @hits = match_positions("CTCTACCCCAACCCC", \$inseq->seq);
print STDERR "found ".scalar(@hits)." hits\n";
foreach (@hits) {
    my ($beg, $end, $hit) = @$_;
    print STDERR "beg: $_->[0]   end: $_->[1]   hit: $_->[2]\n";
}

