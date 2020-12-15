## Reference-based PCR duplicate removal tool

*Cameron Watson*

PCR duplicates are a common issue with Illumina sequencing, since amplification is a necessary step
in Illumina library preparation. These duplicates can cause problems in downstream processes and data analysis.
For instance, having PCR duplicates present can lower the quality of a genome assembly by causing 
convolutions in the DeBruijn graph. Additionally, PCR tends to differentially amplify DNA with different molecular properties. For these reasons, PCR duplicates should be removed from a set of reads; however, this can be computationally intensive if done during early steps of sequence data processing. It is far less computationally intensive to deal with PCR duplicates after alignment, when the data is in Sequence Alginment/Map (SAM) format. 
**This repository contains watson_deduper.py, a tool for removing putative PCR duplicates from a sorted SAM file** 
**of uniquely mapped reads.**

watson_deduper.py can be called from the command-line with the following arguments:

```
--file          A SAM file that has been sorted by RNAME then POS to be deduplicated
--umi           File containing UMIs, unset if randomers
--keepDupes     Specify to write duplicates out to separate sam file, 
                otherwise duplicates are just removed
--correctUMI    Specify to error-correct UMIs that are one nucleotide away from a known UMI; 
                otherwise, reads with unknown UMIs are discarded. When specified, reads with 
                UMIs that are more than one nucleotide away from a known UMI are discarded. 
                Corrected UMIs will appear in the output sam file with an asterisc next to them. 
                Randomers are not error-corrected, and specifying correctUMI without an umi 
                file will result in error.
--help          see full man page for input arguments
```

watson_deduper.py will output a deduplicated SAM file, a PCR duplicates SAM file if keepDupes is specified,
and a summary file containing information including the proportion of duplicate reads and information
regarding UMI filtering and correction if UMI correction was specified.  

Test files and example output files can be found in the test_files directory. 