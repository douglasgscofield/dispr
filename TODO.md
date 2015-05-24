TODO
----

* Allow tags to tag individual primers, and include indication of forward or reverse orientation
* Allow threads in all circumstances
* Indicate primer hits that overlap, perhaps add `--report-overlap`
* Pass full primer hash to match routines?
* multiple primer pairs
* fuzzy matching?
* Begin writing

Thoughts re: addressing multiple issues above:

* how to implement separate tags
  * check for unique tags across forward primers
  * note primer pairs sharing a tag
* how to implement multiplexed primers
  * build a queue of forward primers to search and build a queue of reverse primers to search
  * perhaps just build a queue of primers to search, with orientations part of the queue, this would require making orientation abstracted in the primer info queue
* how to implement flexible use of threads
  * work through a queue of primer searches
* how to implement --orientation
  * do i want this?  i think so... it is easier for the user than post-filtering, it could easily be implemented with a queue... primer orientations are added to the queue, or not


DONE
----

* The `--primers` option may be used to load primer sequences from a Fasta file
* The `--pf` and `--pr` options may be used to specify primer sequences on the command line
* Added `--allow-overlap` option to allow primer hits to overlap
* Assure that searches for hits restart at the proper location when using `--optimise`
* Implemented `--show-mismatches` to count mismatches when `--mismatch-simple` is in effect
* Searches can be limited to specific regions of the reference genome with the `--focal-sites` and `--focal-bounds` options
* Found a way to nicely check for module availability (the whole `re::engine::RE2` issue), need to end the eval with 1: eval 'use re::engine::RE2; 1';
* using re::engine::RE2 for much faster regex matching
* simply reverse-complement the calculated pattern
* Primer pairs delineating each amplicon are indicated.  An amplicon marked `F,R` is produced from the `+` strand, with the forward and reverse primers in native orientation.  An amplicon marked `r,f` is produced from the `-` strand, with both forward and reverse primers in reverse-complement orientation.  An amplicon marked `F,f` was produced from forward primers in native and reverse-complement orientations, and `r,R` was produced from reverse primers.
* Primer hits may be output with `--primer-seq` and `--primer-bed`
* Amplicon hits may be output with `--seq` and `--bed`
* Interior (amplicon minus primers) hits may be output with `--interior-seq` and `--interior-bed`
* The number of unique concrete primer sequences produced with each degenerate primer, including mismatches if specified, is calculated prior to finding matches.  This procedure may be lengthy and can be suppressed with `--skip-count`.
* The `--mismatch-simple` option is provided for specifying mismatches in the 5' end of the primer and optionally in the 3' remainder

