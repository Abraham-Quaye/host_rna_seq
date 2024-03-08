#!/usr/bin/env zsh

filedir=raw_files/annotations
genomedir=raw_files/genome_file
idxdir=raw_files/index

################# COMMAND TO BUILD THEV GENOMIC INDEX #############
echo "Building TURKEY genomic index ..."
hisat2-build -p 8 --ss ${filedir}/turkey_genome.ss --exon ${filedir}/turkey_genome.exons \
${genomedir}/turkey_genome.fa ${idxdir}/turkey_tran &&

echo "Index built successfully" || echo "Program Aborted!!!"
