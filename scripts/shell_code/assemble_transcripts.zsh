#!/usr/bin/env zsh

mapdir=results/hisat2
GTF=raw_files/annotations/turkey_genome_trxpts.gtf
assembled=results/stringtie



for bam in $(ls $mapdir/*.bam)
do
    temp=$(awk -F'.' '{print $1}' <<< "$bam")
    name=$(awk -F'/' '{print $3}' <<< "$temp")
    lab=$(awk -F'hrs' '{print $1}' <<< "${name:7}")
    
    echo "Assembling transcripts for ${name}.bam ..."
    stringtie -p 8 $bam -G $GTF -l $lab -o $assembled/${name:7}.gtf
done


secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Transcript assembly Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs
