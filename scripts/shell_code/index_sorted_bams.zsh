#!/usr/bin/env zsh

mapdir=results/hisat2


for bam in $(ls $mapdir)
do
    echo "Indexing $bam ..."
    samtools index $mapdir/$bam
done

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Read Mapping Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs




