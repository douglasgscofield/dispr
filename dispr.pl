#!/usr/bin/env perl

# TODO: multiple primers
# DONE: remove duplicate hits
# TODO: allow for settable number of mismatches... how?
# TODO: fuzzy matching?
# TODO: identify intervals
# DONE: produce hits for primers only

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
my $o_tag;
my $o_primerhitsonly;
my $o_multiplex = 0;
my $o_orientation = "FR";
my $o_min = 1;
my $o_max = 4000;
my $o_dir = "both";
my $o_ref;
my $o_bed;
my $o_seq;
my $o_verbose;
my $o_help;


my $usage = "
$0  [ OPTIONS ] --ref fasta-file.fa

Applies forward and reverse primer pairs to the reference file,
identifying amplicons within given dimensions.  The primer sequences
are sought paired in specified orientations within a range given
by --min and --max.  The amplicon is measured from the outer extent
of each primer sequence, so includes the primer sequence length.

Primer and search parameters:

    --pf FORWARD_PRIMER   Forward primer sequence (may be specified 2+ times)
    --pr REVERSE_PRIMER   Reverse primer sequence (may be specified 2+ times)
    --pi INT              Index of preloaded primer [default $o_pi]
    --primers FILE        File containing primer pair(s), as forward/reverse
    --tag TAG             String added as tag to output
    --orientation FR      Orientation of primers, only FR supported for now
    --both                Orientation of the reference sequence to search
    --forward             CURRENTLY ONLY --both IS SUPPORTED
    --reverse

Primer hits only:

    --primer-hits-only    Search for primer sequences and write the hit 
                          positions to the bed file and their sequences to the
                          Fasta file.

Amplicons:

    --multiplex           If more than one primer pair presented, consider 
                          amplicons produced by any possible primer pair 
                          (DEFAULT)
    --no-multiplex        If more than one primer pair presented, only consider
                          amplicons produced by each primer pair
    --min bp              Minimum accepted size of amplicon [default $o_min]
    --max bp              Maximum accepted size of amplicon [default $o_max]

Input and output files:

    --ref INPUT_FASTA     Input, Fasta reference sequence in which to find
                          amplicons

    --bed OUTPUT_BED      Output, BED file containing identified amplicon
                          positions
    --seq OUTPUT_FASTA    Output, Fasta sequences containing identified
                          amplicon sequences

Misc:

    --verbose             Describe actions
    --help                Produce this help

";

GetOptions("pf=s"          => \@o_pf,
           "pr=s"          => \@o_pr,
           "pi=i"          => \$o_pi,
           "primers=s"     => \$o_primers,
           "tag=s"         => \$o_tag,
           "orientation=s" => \$o_orientation,
           "both"          => sub { $o_dir = "both" },
           "forward"       => sub { $o_dir = "forward" },
           "reverse"       => sub { $o_dir = "reverse" },
           "primer-hits-only" => \$o_primerhitsonly,
           "multiplex"     => \$o_multiplex,
           "no-multiplex"  => sub { $o_multiplex = 0 },
           "min=i"         => \$o_min,
           "max=i"         => \$o_max,
           "ref=s"         => \$o_ref,
           "bed=s"         => \$o_bed,
           "seq=s"         => \$o_seq,
           "verbose"       => \$o_verbose,
           "help|?"        => \$o_help) or die $usage;
die $usage if $o_help;
#die "only one primer pair currently supported" if @o_pf > 1 or @o_pr > 1;
die "only FR orientation currently supported" if $o_orientation ne "FR";
die "only both strands currently supported" if $o_dir ne "both";
die "only primer hits currently supported" if not $o_primerhitsonly;


sub expand_dot($);  # expand '.' in DNA regex
sub prepare_primer($);  # prepare primer for searches
sub match_positions($$);  # search for primer
sub remove_duplicate_intervals($);  # remove intervals with duplicate beg, end


my %forward = prepare_primer($o_pf[$o_pi]);
my %reverse = prepare_primer($o_pr[$o_pi]);

print STDERR "Assuming primer orientation '$o_orientation' as so, for example primers:

    Forward:F:ACGTCT
    Reverse:R:TTACGC

        Forward>
    ----ACGTCT--------------GCGTAA-------
    ----TGCAGA--------------CGCATT-------
                          <esreveR

Amplicons are identified by being delimited by Forward-esreveR primer pairs,
one in forward orientation, the other in reverse-complement orientation.  Note
that together with their reverse-complements, these primers can delimit
amplicons three additional ways: Reverse-drawroF, Forward-draworF and
Reverse-esreveR.  All of these possibilities are considered here.

In future I may add an indication of which primer pair delimited which
in-silico amplicon.

";

sub iftag  { return $o_tag ? $o_tag    : ""; }
sub iftags { return $o_tag ? "$o_tag " : " "; }

print STDERR iftags()."forward   : $forward{forwardpattern}, $forward{count} expanded sequences\n";
print STDERR iftags()."forward rc: $forward{revcomppattern}, same number in reverse complement\n";
print STDERR iftags()."reverse   : $reverse{forwardpattern}, $reverse{count} expanded sequences\n";
print STDERR iftags()."reverse rc: $reverse{revcomppattern}, same number in reverse complement\n";
print STDERR "\n";

print STDERR "\nProducing output for primer hits only (--primer-hits-only)\n\n" if $o_primerhitsonly;

print STDERR "WARNING: no Fasta or BED output will be produced.\n" if not $o_bed and not $o_seq;

# open input, create output
my ($in, $out_fasta, $out_bed);

$in = Bio::SeqIO->new(-file => "<$o_ref", -format => 'fasta') or 
    die "could not open input file '$o_ref': $!";
if ($o_seq) {
    $out_fasta = Bio::SeqIO->new(-file => ">$o_seq", -format => 'fasta') or
        die "could not open output Fasta file '$o_seq': $!";
}
if ($o_bed) {
    open($out_bed, ">$o_bed") or 
        die "could not open output BED file '$o_bed': $!";
}

my $Un;

while (my $inseq = $in->next_seq()) {
    my $seqname = $inseq->display_id();
    if ($seqname =~ /^chrUn/) {
        print STDERR "chrUn*:".iftag()." searching ...\n" if not $Un;
        ++$Un;
    } else {
        print STDERR "$seqname:".iftag()." searching ...\n";
    }
    # any _forward_hits can be complemented by any _revcomp_hits
    my @f_forward_hits = match_positions($forward{forwardquoted}, \$inseq->seq);
    my @r_revcomp_hits = match_positions($reverse{revcompquoted}, \$inseq->seq);
    my @r_forward_hits = match_positions($reverse{forwardquoted}, \$inseq->seq);
    my @f_revcomp_hits = match_positions($forward{revcompquoted}, \$inseq->seq);

    # sort and remove duplicate hits that start at the same position
    my @forward_hits = sort { $a->[0] <=> $b->[0] } (@f_forward_hits, @r_forward_hits);
    my @revcomp_hits = sort { $a->[0] <=> $b->[0] } (@f_revcomp_hits, @r_revcomp_hits);
    my $n_forward_dups = remove_duplicate_intervals(\@forward_hits);
    my $n_revcomp_dups = remove_duplicate_intervals(\@revcomp_hits);

    next if not @forward_hits and not @revcomp_hits;  # no hits found

    print STDERR $seqname.":".iftag().
                 " forward hits ".scalar(@forward_hits)." ($n_forward_dups dups),".
                 " revcomp hits ".scalar(@revcomp_hits)." ($n_revcomp_dups dups)\n";

    if ($o_primerhitsonly) {
        # join hits, sort, produce bed intervals and Fasta file
        my @all_hits = sort { $a->[0] <=> $b->[0] } ( @forward_hits, @revcomp_hits );
        my $n_all_dups = remove_duplicate_intervals(\@all_hits);
        print STDERR "$seqname:".iftag()." all hits ".scalar(@all_hits)." ($n_all_dups dups)\n";

        foreach my $h (@all_hits) {
            if ($o_bed) {
                my $name = $h->[2];  # the hit sequence itself
                $name = "$o_tag:$name" if $o_tag;
                print $out_bed $seqname."\t".$h->[0]."\t".$h->[1]."\t".$name."\n";
            }
            if ($o_seq) {
                # use base-1 GFF-type intervals in Fasta name
                my $name = "$seqname:".($h->[0] + 1)."-".$h->[1];
                $name .= ":$o_tag" if $o_tag;
                my $hitseq = Bio::Seq->new(-id => $name,
                                           -seq => $h->[2],
                                           -alphabed => 'dna');
                $out_fasta->write_seq($hitseq);
            }
        }

        next;  # primers only, skip to the next input sequence
    }

    # calculate amplicon intervals, produce bed intervals
    #foreach (@all_hits) {
    #    my ($name, $tag, $beg, $end, $hit) = @$_;
    #    print STDERR "main: $name   $tag   $beg-$end   $hit\n";
    #}
}

$in->close();
$out_fasta->close();
$out_bed->close();


# ---- local subroutines ----------------------------------- 



# By default Bio::Tools::SeqPattern->expand() replaces N with . in the
# regex it returns.  Though this is strictly correct when enforcing a
# DNA alphabet, I would prefer it be more explicit.  This replaces '.'
# with '[ACTGN]'.
#
sub expand_dot($) {
    my $pat = shift;
    $pat =~ s/\./[ACTGN]/g;
    return $pat;
}



# Prepare a primer for searching by creating a hash for a primer sequence,
# passed as the single argument in one of three forms:
#
#    1. straight sequence:    CTYNARG...
#    2. annotated sequence:   name:direction:CTYNARG...
#    3. Bio::Seq object (NOT IMPLEMENTED YET)
#
# it prepares a hash containing objects prepared for searching for that primer.
# Included are:
#
#     name            name of the primer (if supplied)
#     dir             direction of the primer (if supplied)
#     sequence        sequence of the primer as provided
#     Seq             Bio::Seq object for the primer
#     SeqPattern      Bio::Tools::SeqPattern object for the primer
#     forwardpattern  regex for the forward (given) orientation of the primer
#     forwardquoted   a quoted 'qr/.../i' version of forwardpattern
#     revcomppattern  regex for the reverse complement of the given primer
#     revcompquoted   a quoted 'qr/.../i' version of revcomppattern
#     IUPAC           Bio::Tools::IUPAC object for the given primer
#     count           number of concrete sequences formable from primer
#
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



# Passed in a pattern quoted with 'qr/.../i' and a reference to a sequence to
# search.  Returns an array of anonymous arrays containing the 0-based
# beginning and end of the hit and the sequence of the hit.  The interval is
# [beg, end), the same as a BED interval, and each anonymous array contains
#
# [ $beg, $end, $hit_sequence ]
#
sub match_positions($$) {
    my ($pat, $seq) = @_;
    my @ans;
    while ($$seq =~ /$pat/ig) {
        my ($beg, $end) = ($-[0], $+[0]);
        my $hit = substr($$seq, $beg, $end - $beg);
        print STDERR "match_positions: $beg-$end   $hit\n" if $o_verbose;
        push @ans, [ $beg, $end, $hit ];
    }
    return @ans;
}



# When pass a reference to an array of intervals, removes all intervals with the
# same start and stop sites (->[0] and ->[1]) and returns the number of duplicate
# intervals removed.
#
sub remove_duplicate_intervals($) {
    my $a = shift;
    my %seen;
    my $n = scalar(@$a);
    @$a = grep { ! $seen{$_->[0]."-".$_->[1]}++ } @$a;
    return $n - scalar(@$a);
}



