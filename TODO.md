TODO
----

* find a way to nicely check for module availability (the whole `re::engine::RE2` issue)
* multiple primer pairs
* primer pairs from Fasta file
* fuzzy matching?


DONE
----

* using re::engine::RE2 for much faster regex matching
* simply reverse-complement the calculated pattern
* Identify amplicons by the primer pair making them... so + strand and - strand are distinguishable
* --interior-bed and --interior-seq to output amplicons minus primers
* count the number of concrete primer sequences with --mismatch-simple
* settable number of mismatches at 5' end with --mismatch-simple
* remove duplicate hits
* identify amplicon intervals
* produce hits for primers only

