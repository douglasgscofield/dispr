Degenerate In-Silico PcR
========================

Given a pair of degenerate primers, find sites in the given genome matching some concrete primer and produce one or all of six types of output files:

* Fasta-format file of primer-pair amplicon sequences
* BED-format file of primer-pair amplicons
* Fasta-format file of primer binding site sequences
* BED-format file of primer binding sites
* Fasta-format file of interior sequences, amplicons minus primers
* BED-format file of interior sequences, amplicons minus primers

Some mismatch sites can be allowed within the primer.  Following typical primer binding dynamics, the 5' "head" of the primer can have up to 5 mismatches over a user-specified length, with the remaining 3' "tail" of the primer holding up to 2 mismatches.  The number of concrete sequences arising from degenerate primers is also calculated, though if a number of mismatches are allowed, this can take several minutes because an explicit counting method is used to accomodate redundancy arising from allowing for mismatches in already degenerate sequences; the count can be disabled with the `--skip-count` option.

Dispr is implemented using Perl regular expression matching.  The `re::engine::RE2` regular expression drop-in module is used if available as it is much faster than Perl's implementation for the types of patterns matched here.  If you are having trouble compiling `re::engine::RE2` especially on Mac OS X, try following the recommendations [here](http://stackoverflow.com/questions/17440137/error-variable-length-array-of-non-pod-element-type-when-compiling-reengine/29755725#29755725).



Options
-------

Full help is available with `dispr --help`.  Some options of note:

`--mismatch-simple INT1:INT2[:INT3]`
: Allow for *INT1* mismatches in the 5'-most *INT2* bases ("head") of each primer, with optionally *INT3* mismatches in the remaining 3' portion ("tail") of the primer

`--show-mismatches`
: Count the number of mismatches per primer hit when `--mismatch-simple` is in effect

`--optimise`
: Try to speed up searches with `--mismatch-simple` by searching for matches against the (presumably less-complex) tail portion of each primer first, and only then looking for a match against the adjacent primer head

`--overlap`
: Allow primer hits to overlap; amplicons may overlap at any time if produced by nonoverlapping primer hits

`--focal-sites BED`
: Confine primer searches to regions of the reference genome described in the BED-format file.  The option `--focal-bounds` allows the specification of a further region up- and downstream of each BED interval that is included in the search (default &plusmn;1000 bp)

`--threads N`
: Execute *N* (2 or 4) primer searches in parallel


License
-------
Gnu Public License v2

