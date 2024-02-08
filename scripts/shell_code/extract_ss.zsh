#!/usr/bin/env zsh

filedir=raw_files/annotations

################# EXTRACT SPLICESITES FROM THEV GTF FILE ################
echo "extracting TURKEY splice-sites ..."
hisat2_extract_splice_sites.py $filedir/turkey_genome.gtf > $filedir/turkey_genome.ss &&
echo "TURKEY splice-site extraction completed successfully"

