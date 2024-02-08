#!/usr/bin/env zsh

mapdir=results/hisat2
GTF=raw_files/annotations/turkey_genome.gtf
assembled=results/stringtie



for bam in $(ls $mapdir)
do
    name=$(awk -F'.' '{print $1}' <<< "$bam")
    lab=$(awk -F'hrs' '{print $1}' <<< "${name:7}")
    
    echo "Assembling transcripts for $bam ..."
    stringtie -p 8 $mapdir/$bam -G $GTF -l $lab -o $assembled/${name:7}.gtf
done


secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Read Mapping Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs
