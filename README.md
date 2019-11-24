# find_circ_any
Tools to identify circRNAs originating from any genomic loci from RNA-seq data

find_circ_any.py is used to identify circRNAs from single RNA-seq.

find_circ_pair.py is used to identify circRNAs from paired-end RNA-seq.

Usage:
python find_circ_any.py candidate.sam sites.bed sites.reads -G /path_to_your_chromosome_folder -p prefix -s sites.log
