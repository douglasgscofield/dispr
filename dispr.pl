#!/usr/bin/env perl

# TODO: multiple primers
# TODO: remove duplicate hits
# TODO: allow for settable number of mismatches... how?
# TODO: fuzzy matching?
# TODO: identify intervals

use strict;
use warnings;

use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;
use Bio::Tools::SeqPattern;

#"class1:F:CTNCAYVARCCYATGTAYYWYTTBYT"
#"class1:R:GTYYTNACHCYRTAVAYRATRGGRTT"

#my @o_pf = qw/ class1:F:CTNCAYVARCCYATGTAYYWYTTBYT class1:R:CTNCANWCNCCHATGTAYTTYYTBCT /;
#my @o_pr = qw/ class2:F:GTYYTNACHCYRTAVAYRATRGGRTT class2:R:TTYCTBARRSTRTARATNADRGGRTT /;
my @o_pf = qw/ CTNCAYVARCCYATGTAYYWYTTBYT CTNCANWCNCCHATGTAYTTYYTBCT /;
my @o_pr = qw/ GTYYTNACHCYRTAVAYRATRGGRTT TTYCTBARRSTRTARATNADRGGRTT /;
my $o_pi = 0;
#>class1|F
#CTNCAYVARCCYATGTAYYWYTTBYT
#>class1|R
#GTYYTNACHCYRTAVAYRATRGGRTT
#>class2|F
#CTNCANWCNCCHATGTAYTTYYTBCT
#>class2|R
#TTYCTBARRSTRTARATNADRGGRTT
my $o_primers;
my $o_multiplex = 0;
my $o_orientation = "FR";
my $o_min = 1;
my $o_max = 4000;
my $o_dir = "both";
my $o_ref = "test.fa";
my $o_bed = "out.bed";
my $o_seq = "out.fa";
my $o_verbose;
my $o_help;


my $usage = "
$0  --forward FORWARD_PRIMER  --reverse REVERSE_PRIMER --orientation FR|RF|FF|ALL  --reference FASTA_FILE

Applies forward and reverse primer pairs to the reference file,
identifying amplicons within given dimensions.  The primer sequences
are sought paired in specified orientations within a range given
by --min and --max.  The amplicon is measured from the outer extent
of each primer sequence, so includes the primer sequence length.

    --pf FORWARD_PRIMER   Forward primer sequence (may be specified more than once)
    --pr REVERSE_PRIMER   Reverse primer sequence (may be specified more than once)
    --pi INT              Index of preloaded primer [default $o_pi]
    --primers FILE        File containing primer pair(s), as forward/reverse
    --multiplex           If more than one primer pair presented, consider amplicons
                          produced by any possible primer pair (DEFAULT)
    --no-multiplex        If more than one primer pair presented, only consider
                          amplicons produced by each primer pair
    --orientation FR      Orientation of primers, only FR supported for now
    --min bp              Minimum accepted size of amplicon [default $o_min]
    --max bp              Maximum accepted size of amplicon [default $o_max]

    --both                Orientation of the reference sequence
    --forward
    --reverse

    --ref INPUT_FASTA     Input, Fasta reference sequence in which to find amplicons

    --bed OUTPUT_BED      Output, BED file containing identified amplicon positions
    --seq OUTPUT_FASTA    Output, Fasta sequences containing identified amplicon sequences

    --verbose             Describe actions
    --help                Produce this help

";

GetOptions("pf=s"          => \@o_pf,
           "pr=s"          => \@o_pr,
           "pi=i"          => \$o_pi,
           "primers=s"     => \$o_primers,
           "multiplex"     => \$o_multiplex,
           "no-multiplex"  => sub { $o_multiplex = 0 },
           "orientation=s" => \$o_orientation,
           "min=i"         => \$o_min,
           "max=i"         => \$o_max,
           "both"          => sub { $o_dir = "both" },
           "forward"       => sub { $o_dir = "forward" },
           "reverse"       => sub { $o_dir = "reverse" },
           "ref=s"         => \$o_ref,
           "bed=s"         => \$o_bed,
           "seq=s"         => \$o_seq,
           "verbose"       => \$o_verbose,
           "help|?"        => \$o_help) or die $usage;
die $usage if $o_help;
#die "only one primer pair currently supported" if @o_pf > 1 or @o_pr > 1;
die "only FR orientation currently supported" if $o_orientation ne "FR";
die "only both strands currently supported" if $o_dir ne "both";

sub expand_dot($) {
    # would rather be explicit about '.'
    my $pat = shift;
    $pat =~ s/\./[ACTGN]/g;
    return $pat;
}

# http://stackoverflow.com/questions/87380/how-can-i-find-the-location-of-a-regex-match-in-perl
sub match_positions($$) {
    my ($pat, $seq) = @_;
    my @ans;
    while ($$seq =~ /$pat/ig) {
        my ($beg, $end) = ($-[0], $+[0]);
        my $hit = substr($$seq, $beg, $end - $beg);
        print STDERR "match_positions: $beg-$end   $hit\n" if $o_verbose;
        push @ans, [$beg, $end, $hit];
    }
    return @ans;
}

sub prepare_primer($) {
    my ($primer) = @_;
    my %dest;
    if (index($primer, ":") >= 0) {
        my @p = split /:/, $primer;
        die "format is   name:dir:sequence" if @p != 3;
        $dest{name} = $p[0];
        $dest{dir} = $p[1];
        $primer = $p[2];
    }
    $dest{sequence} = $primer;
    my $s = Bio::Seq->new(-seq => $primer, -alphabet => 'dna');
    $dest{Seq} = $s;
    my $seqpattern = Bio::Tools::SeqPattern->new(-seq => $s->seq(), -type => 'dna');
    $dest{SeqPattern} = $seqpattern;
    $dest{forwardpattern} = expand_dot($seqpattern->expand());
    $dest{forwardquoted} = qr/$dest{forwardpattern}/i;
    $dest{revcomppattern} = expand_dot($seqpattern->revcom(1)->expand());
    $dest{revcompquoted} = qr/$dest{revcomppattern}/i;
    my $iupac = Bio::Tools::IUPAC->new(-seq => $s);
    $dest{IUPAC} = $iupac;
    $dest{count} = $iupac->count();
    return %dest;
}


my %forward = prepare_primer($o_pf[$o_pi]);
my %reverse = prepare_primer($o_pr[$o_pi]);

print STDERR "forward: $forward{forwardpattern}, $forward{count} expanded sequences\n";
print STDERR "reverse: $reverse{forwardpattern}, $reverse{count} expanded sequences\n";

my $in = Bio::SeqIO->new(-file => "<$o_ref", -format => 'fasta');

my $Un;

while (my $inseq = $in->next_seq()) {
    my $seqname = $inseq->display_id();
    if ($seqname =~ /^chrUn/) {
        print STDERR "chrUn* searching ...\n" if not $Un;
        ++$Un;
    } else {
        print STDERR "$seqname: searching ...\n";
    }
    # any _forward_hits can be complemented by any _revcomp_hits
    my @f_forward_hits = match_positions($forward{forwardquoted}, \$inseq->seq);
    my @r_revcomp_hits = match_positions($reverse{revcompquoted}, \$inseq->seq);
    my @f_revcomp_hits = match_positions($forward{revcompquoted}, \$inseq->seq);
    my @r_forward_hits = match_positions($reverse{forwardquoted}, \$inseq->seq);

    # sort and remove duplicate hits that start at the same position
    my @forward_hits = sort { $a->[0] <=> $b->[0] } (@f_forward_hits, @r_forward_hits);
    my $n = @forward_hits;
    my %seen;
    @forward_hits = grep { ! $seen{$_->[0]}++ } @forward_hits;
    my $n_forward_dups = $n - @forward_hits;
    my @revcomp_hits = sort { $a->[0] <=> $b->[0] } (@f_revcomp_hits, @r_revcomp_hits);
    $n = @revcomp_hits;
    undef %seen;
    @revcomp_hits = grep { ! $seen{$_->[0]}++ } @revcomp_hits;
    my $n_revcomp_dups = $n - @revcomp_hits;

    #my @hits = match_positions($seqname, "AACCCTACCTAAACCTCA", \$inseq->seq);
    #my @hits = match_positions($seqname, "CTCTACCCCAACCCC", \$inseq->seq);
    if ($o_verbose) {
        print STDERR "$seqname found ".scalar(@forward_hits)." forward hits (".
                     $n_forward_dups, " dups), and ".scalar(@revcomp_hits).
                     " revcomp hits (".$n_revcomp_dups." dups)\n";
    }
    #foreach (@all_hits) {
    #    my ($name, $tag, $beg, $end, $hit) = @$_;
    #    print STDERR "main: $name   $tag   $beg-$end   $hit\n";
    #}
}


