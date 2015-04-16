#!/usr/bin/env perl

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
my $o_pi = 1;  # 1 - 1 = 0
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
my $o_primerbed;
my $o_primerseq;
my $o_multiplex = 0;
my $o_orientation = "FR";
my $o_min = 1;
my $o_max = 1200;
my $o_maxmax = 10000;
my $o_dir = "both";
my $o_ref;
my $o_bed;
my $o_seq;
my $o_expand_dot = 0;
my $o_verbose;
my $o_help;

my $o_mismatch_simple;
my $o_mm_int1;
my $o_mm_int2;


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
         1: class 1 f $o_pf[0]  r $o_pr[0]
         2: class 2 f $o_pf[1]  r $o_pr[1]
    --primers FILE        File containing primer pair(s), as forward/reverse
                          (CURRENTLY UNSUPPORTED)
    --tag TAG             String added as tag to output
    --orientation FR      Orientation of primers, only FR supported for now
    --both                Orientation of the reference sequence to search
    --forward             CURRENTLY ONLY --both IS SUPPORTED
    --reverse

    --mismatch-simple INT1:INT2
                          Allow up to INT1 mismatches in 5'-most INT2 bp of 
                          each primer.  INT1 must be <= 3 and INT2 must be
                          <= 7.  To allow up to 2 mismatches in the 5'-most
                          5 bp:  --mismatch-simple 2:5

Amplicons:

    --multiplex           If more than one primer pair presented, consider 
                          amplicons produced by any possible primer pair 
                          (DEFAULT, BUT >1 PRIMER PAIR NOT CURRENTLY SUPPORTED)
    --no-multiplex        If more than one primer pair presented, only consider
                          amplicons produced by each primer pair
    --min bp              Minimum accepted size of amplicon [default $o_min]
    --max bp              Maximum accepted size of amplicon [default $o_max]
    --maxmax bp           Maximum 'too-long' amplicon to track [default $o_maxmax]
    --output-extremes     Include too-short and too-long amplicons in output,
                          the names include resp. SHORT and LONG
                          (NOT CURRENTLY SUPPORTED)

Input and output files:

    --ref INPUT_FASTA     Input, Fasta reference sequence in which to find
                          amplicons

    --bed OUTPUT_BED      Output, BED file containing identified amplicon
                          positions
    --seq OUTPUT_FASTA    Output, Fasta sequences containing identified
                          amplicon sequences
    --primer-bed BED      Output, BED file containing hits for primer
                          sequences found
    --primer-seq FASTA    Output, Fasta sequences containing identified
                          primer sequences

Misc:

    --expand-dot          Expand '.' in regexs to '[ACGTN]'
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
           "mismatch-simple=s" => \$o_mismatch_simple,
           "primer-bed=s"  => \$o_primerbed,
           "primer-seq=s"  => \$o_primerseq,
           "multiplex"     => \$o_multiplex,
           "no-multiplex"  => sub { $o_multiplex = 0 },
           "min=i"         => \$o_min,
           "max=i"         => \$o_max,
           "maxmax=i"      => \$o_maxmax,
           "ref=s"         => \$o_ref,
           "bed=s"         => \$o_bed,
           "seq=s"         => \$o_seq,
           "verbose"       => \$o_verbose,
           "expand-dot"    => \$o_expand_dot,
           "help|?"        => \$o_help) or die $usage;
die $usage if $o_help;
#die "only one primer pair currently supported" if @o_pf > 1 or @o_pr > 1;
die "only FR orientation currently supported" if $o_orientation ne "FR";
die "only both strands currently supported" if $o_dir ne "both";
die "must provide results name --tag" if not $o_tag;
die "must provide sequence to search with --ref" if not $o_ref;

if ($o_mismatch_simple) {
    ($o_mm_int1, $o_mm_int2) = split(/:/, $o_mismatch_simple, 2);
    die "unable to interpret --mismatch-simple argument" if not $o_mm_int1 or not $o_mm_int2;
    die "must allow <= 3 mismatches" if $o_mm_int1 > 3;
    die "must span <= 7 5' bases" if $o_mm_int2 > 7;
}

sub expand_dot($);  # expand '.' in DNA regex
sub prepare_primer($);  # prepare primer for searches
sub apply_mismatch_simple($$$$);  # prepare query sequence for mismatches
sub match_positions($$);  # search for primer
sub remove_duplicate_intervals($);  # remove intervals with duplicate beg, end
sub dump_primer_hits($$$);  # dump primer-only intervals
sub dump_amplicon_hits($$$$);  # calculate and dump amplicons


my %forward = prepare_primer($o_pf[$o_pi - 1]);
my %reverse = prepare_primer($o_pr[$o_pi - 1]);

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

print STDERR "WARNING: no Fasta or BED output will be produced.\n" if not $o_bed and not $o_seq and not $o_primerbed and not $o_primerseq;

# open input, create output
my ($in, $out_seq, $out_bed, $out_primerseq, $out_primerbed);

$in = Bio::SeqIO->new(-file => "<$o_ref", -format => 'fasta') or 
    die "could not open input file '$o_ref': $!";
if ($o_seq) {
    $out_seq = Bio::SeqIO->new(-file => ">$o_seq", -format => 'fasta') or
        die "could not open output Fasta file '$o_seq': $!";
}
if ($o_bed) {
    open($out_bed, ">$o_bed") or 
        die "could not open output BED file '$o_bed': $!";
}
if ($o_primerseq) {
    $out_primerseq = Bio::SeqIO->new(-file => ">$o_primerseq", -format => 'fasta') or
        die "could not open output Fasta file '$o_primerseq': $!";
}
if ($o_primerbed) {
    open($out_primerbed, ">$o_primerbed") or 
        die "could not open output BED file '$o_primerbed': $!";
}

my $Un;

while (my $inseq = $in->next_seq()) {
    my $this_seqname = $inseq->display_id();
    my $this_sequence = $inseq->seq();
    if ($this_seqname =~ /^chrUn/) {
        print STDERR "chrUn*:".iftag()." searching ...\n" if not $Un;
        ++$Un;
    } else {
        print STDERR "$this_seqname:".iftag()." searching ...\n";
    }
    # any _forward_hits can be complemented by any _revcomp_hits
    my @f_forward_hits = match_positions($forward{forwardquoted}, \$this_sequence);
    my @r_revcomp_hits = match_positions($reverse{revcompquoted}, \$this_sequence);
    my @r_forward_hits = match_positions($reverse{forwardquoted}, \$this_sequence);
    my @f_revcomp_hits = match_positions($forward{revcompquoted}, \$this_sequence);

    # sort and remove duplicate hits that start at the same position
    my @forward_hits = sort { $a->[0] <=> $b->[0] } (@f_forward_hits, @r_forward_hits);
    my @revcomp_hits = sort { $a->[0] <=> $b->[0] } (@f_revcomp_hits, @r_revcomp_hits);
    my $n_forward_dups = remove_duplicate_intervals(\@forward_hits);
    my $n_revcomp_dups = remove_duplicate_intervals(\@revcomp_hits);

    next if not @forward_hits and not @revcomp_hits;  # no hits found

    print STDERR $this_seqname.":".iftag().
                 " forward hits ".scalar(@forward_hits)." ($n_forward_dups dups),".
                 " revcomp hits ".scalar(@revcomp_hits)." ($n_revcomp_dups dups)\n";

    dump_primer_hits($this_seqname, \@forward_hits, \@revcomp_hits);

    dump_amplicon_hits($this_seqname, \$this_sequence, \@forward_hits, \@revcomp_hits);

}

$in->close() if $in;
$out_seq->close() if $out_seq;
$out_bed->close() if $out_bed;
$out_primerseq->close() if $out_primerseq;
$out_primerbed->close() if $out_primerbed;


# ---- local subroutines ----------------------------------- 



# By default Bio::Tools::SeqPattern->expand() replaces N with . in the
# regex it returns.  Though this is strictly correct when enforcing a
# DNA alphabet, I would prefer it be more explicit.  This replaces '.'
# with '[ACTGN]'.
#
sub expand_dot($) {
    my $pat = shift;
    $pat =~ s/\./[ACTGN]/g if $o_expand_dot;
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
    if ($o_mismatch_simple) {
        $dest{forwardpattern} = apply_mismatch_simple($primer, $o_mm_int1, $o_mm_int2, 0);
        $dest{revcomppattern} = apply_mismatch_simple($primer, $o_mm_int1, $o_mm_int2, 1);
        $dest{mismatch} = $o_mismatch_simple;
        $dest{forward0} = expand_dot($seqpattern->expand());
        $dest{revcomp0} = expand_dot($seqpattern->revcom(1)->expand());
        print STDERR "forwardpattern    $dest{forwardpattern}\n" if $o_verbose;
        print STDERR "revcomppattern    $dest{revcomppattern}\n" if $o_verbose;
        print STDERR "forward0          $dest{forward0}\n" if $o_verbose;
        print STDERR "revcomp0          $dest{revcomp0}\n" if $o_verbose;
    } else {
        $dest{forwardpattern} = expand_dot($seqpattern->expand());
        $dest{revcomppattern} = expand_dot($seqpattern->revcom(1)->expand());
    }
    $dest{forwardquoted} = qr/$dest{forwardpattern}/i;
    $dest{revcompquoted} = qr/$dest{revcomppattern}/i;
    my $iupac = Bio::Tools::IUPAC->new(-seq => $s);
    $dest{IUPAC} = $iupac;
    $dest{count} = $iupac->count();
    return %dest;
}



sub apply_mismatch_simple($$$$) {
    my ($p, $m, $len, $is_rc) = @_;
    print STDERR "p = $p, m = $m, len = $len, is_rc = $is_rc\n";
    my ($head, $tail);
    #if ($is_rc) {
    #    $head = substr($p, -$len);   $tail = substr($p, 0, length($p) - $len);
    #} else {
        $head = substr($p, 0, $len); $tail = substr($p, $len);
    #}
    die "no tail available" if not $tail;
    print STDERR "head = $head, tail = $tail\n" if $o_verbose;
    my @p;
    if ($m == 3) {
        for (my $i = 0; $i < $len - 2; ++$i) {
            my $s = $head;
            substr($s, $i, 1) = "N";
            for (my $j = $i + 1; $j < $len - 1; ++$j) {
                my $ss = $s;
                substr($ss, $j, 1) = "N";
                for (my $k = $j + 1; $k < $len; ++$k) {
                    my $sss = $ss;
                    substr($sss, $k, 1) = "N";
                    my $sp = Bio::Tools::SeqPattern->new(-seq => $sss, -type => 'dna');
                    my $pat = expand_dot($is_rc ? $sp->revcom(1)->expand() : $sp->expand());
                    print STDERR "i = $i, pat = $pat\n" if $o_verbose;
                    push @p, $pat;
                }
            }
        }
    } elsif ($m == 2) {
        for (my $i = 0; $i < $len - 1; ++$i) {
            my $s = $head;
            substr($s, $i, 1) = "N";
            for (my $j = $i + 1; $j < $len; ++$j) {
                my $ss = $s;
                substr($ss, $j, 1) = "N";
                my $sp = Bio::Tools::SeqPattern->new(-seq => $ss, -type => 'dna');
                my $pat = expand_dot($is_rc ? $sp->revcom(1)->expand() : $sp->expand());
                print STDERR "i = $i, pat = $pat\n" if $o_verbose;
                push @p, $pat;
            }
        }
    } elsif ($m == 1) {
        for (my $i = 0; $i < $len; ++$i) {
            my $s = $head;
            substr($s, $i, 1) = "N";
            my $sp = Bio::Tools::SeqPattern->new(-seq => $s, -type => 'dna');
            my $pat = expand_dot($is_rc ? $sp->revcom(1)->expand() : $sp->expand());
            print STDERR "i = $i, pat = $pat\n" if $o_verbose;
            push @p, $pat;
        }
    } else {
        die "unrecognised number of mismatches $m";
    }
    @p = reverse @p if $is_rc;
    my $m_pat = '(' . join('|', @p) . ')';
    my $tp = Bio::Tools::SeqPattern->new(-seq => $tail, -type => 'dna');
    my ($tailpat, $fullpat);
    if ($is_rc) {
        $tailpat = expand_dot($tp->revcom(1)->expand());
        $fullpat = $tailpat . $m_pat;
    } else {
        $tailpat = expand_dot($tp->expand());
        $fullpat = $m_pat . $tailpat;
    }
    print STDERR "m_pat = $m_pat, tailpat = $tailpat, fullpat = $fullpat\n" if $o_verbose;
    return $fullpat;
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



# Dump hits for primer sequences, if requested.  Join hits, sort,
# produce BED and/or Fasta files.
#
sub dump_primer_hits($$$) {
    my ($seqname, $forward_hits, $revcomp_hits) = @_;

    my @all_hits = sort { $a->[0] <=> $b->[0] } ( @$forward_hits, @$revcomp_hits );
    my $n_all_dups = remove_duplicate_intervals(\@all_hits);
    print STDERR "dump_primer_hits: $seqname:".iftag()." all hits ".
                 scalar(@all_hits)." ($n_all_dups dups)\n" if $o_verbose;

    return if not $o_primerbed and not $o_primerseq;

    foreach my $h (@all_hits) {
        if ($o_primerbed) {
            my $name = $h->[2];  # the hit sequence itself
            $name = "$o_tag:$name" if $o_tag;
            #print $out_primerbed $seqname."\t".$h->[0]."\t".$h->[1]."\t".$name."\n";
            $out_primerbed->print($seqname."\t".$h->[0]."\t".$h->[1]."\t".$name."\n");
        }
        if ($o_primerseq) {
            # use base-1 GFF-type intervals in Fasta name
            my $name = "$seqname:".($h->[0] + 1)."-".$h->[1];
            $name .= ":$o_tag" if $o_tag;
            my $hitseq = Bio::Seq->new(-id => $name,
                                       -seq => $h->[2],
                                       -alphabed => 'dna');
            $out_primerseq->write_seq($hitseq);
        }
    }
}



# Dump hits for amplicon sequences, if requested.  Construct amplicons,
# sort into too-short, too-long, and just-right lengths, count them up,
# and produce BED and/or Fasta files.
#
sub dump_amplicon_hits($$$$) {
    my ($seqname, $sequence, $forward_hits, $revcomp_hits) = @_;
    # intervals extend between forward_hits->[0] and revcomp_hits->[1]
    # calculate amplicon intervals, produce bed intervals
    my @amp_short; # too short  0 <=      < $o_min
    my @amp_long; # too long   $o_max <  <= $o_maxmax
    my @amp; # just right $o_min <= <= $o_max
    foreach my $f (@$forward_hits) {
        foreach my $r (@$revcomp_hits) {
            next if $r->[0] < $f->[1];  # at least de-overlap the primers
            my ($beg, $end) = ( $f->[0], $r->[1] );
            my $amp_len = $end - $beg;
            next if $amp_len < 0 or $amp_len > $o_maxmax;
            my $amplicon = substr($$sequence, $beg, $end - $beg);
            my $arr = [ $beg, $end, $amplicon ];
            if ($amp_len < $o_min) {
                push @amp_short, $arr;
            } elsif ($amp_len > $o_max) {
                push @amp_long, $arr;
            } else {
                push @amp, $arr;
            }
        }
    }
    my $n_dup = remove_duplicate_intervals(\@amp);
    my $n_short_dup = remove_duplicate_intervals(\@amp_short);
    my $n_long_dup = remove_duplicate_intervals(\@amp_long);
    print STDERR $seqname.":".iftag().
                 " right ".scalar(@amp)." ($n_dup dups),".
                 " tooshort ".scalar(@amp_short)." ($n_short_dup dups),".
                 " toolong ".scalar(@amp_long)." ($n_long_dup dups)\n";

    return if not $o_bed and not $o_seq;

    foreach my $h (@amp) {
        if ($o_bed) {
            my $name = "$o_tag:".length($h->[2]);
            #print $out_bed $seqname."\t".$h->[0]."\t".$h->[1]."\t".$name."\n";
            $out_bed->print($seqname."\t".$h->[0]."\t".$h->[1]."\t".$name."\n");
        }
        if ($o_seq) {
            # use base-1 GFF-type intervals in Fasta name
            my $name = "$seqname:".($h->[0] + 1)."-".$h->[1].":$o_tag:".length($h->[2]);
            my $hitseq = Bio::Seq->new(-id => $name,
                                       -seq => $h->[2],
                                       -alphabed => 'dna');
            $out_seq->write_seq($hitseq);
        }
    }
}
