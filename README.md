Degenerate In-Silico PcR
========================

Given a pair of degenerate primers, find sites in the given genome matching some concrete primer and produce one or all of six types of output files:

* Fasta-format file of primer-pair amplicon sequences
* BED-format file of primer-pair amplicons
* Fasta-format file of primer binding site sequences
* BED-format file of primer binding sites
* Fasta-format file of interior sequences, amplicons minus primers
* BED-format file of interior sequences, amplicons minus primers

A small number of mismatches can be allowed at the 5' end of the primer.  The number of concrete sequences arising from degenerate primers, including when mismatches are allowed, is also calculated.

License
-------
Gnu Public License v2
