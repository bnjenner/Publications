
# HTSeq-TAG-Counts

HTSeq is a set of tools written by Simon Anders (https://github.com/simon-anders/htseq) for processing and analyzing data from high-throughput sequencing experiments. This modified version of the "htseq-counts" function was written with quantifying TAGseq reads for non-model orgranisms in mind. Keep in mind, this program is 95% Simon Anders, these modifications were necessary for TAGseq experiments conducted in the Gordon Lab by Bradley N. Jenner and Peter M. Henry.  

## Purpose
Annotation files for non-model organisms can often be unreliable or incomplete. For this reason, manual annotations are often necessary. However, this can potentially generate large numbers of overlapping features in the annotation files. The default version of "htseq-counts" marks reads mapping to these regions as ambiguous, leading to inaccurate read quantification. In a TAGseq situation where many of the features to be quantified are CDS regions and manually added 3' UTRs, the potential for these features to overlap with downstream genes is high. `tag_counts.py` implements several features that resolve this issue with TAGseq specific issue. 

## Modifications 
Based on the HTSeq (version=0.11.2) "counts.py" program (python3), this script has been modified for more accurate quantification of these potentially overlapping TAGseq reads in non-model orgranisms. This includes the addition of a "--utr-tag" parameter which allows you to specify UTR specific feature annotations that the program will count in addition to specified features (ie. CDS). In situations where reads map to two overlaping features, instead of being marked as ambiguous, "tag_counts.py" will check to see if one of the overlaps is marked with the annotation specifed by the "--utr-tag" option. If this is the case, the read will be counted towards the UTR marked feature, otherwise it will be marked as ambiguous. This modification significantly reduces the number of reads marked ambiguous by "htseq-counts," many of which were generated artifically by manual annotating the GFF/GTF files. 

## Installation 
This script requires that the actual HTSeq program (version 0.11.2) be installed on your computer. Version specific installation can be found at the original repository (https://github.com/simon-anders/htseq), but installation with conda is recommended.

  `conda install -c bioconda HTSeq=0.11.2`
  
To install this program, clone this repository.

## Usage 
  `./tag_counts.py --utr-tag="YOUR_UTR_TAG" INPUT.BAM INPUT.GTF`
