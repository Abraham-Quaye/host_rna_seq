#!/usr/bin/env bash

mapdir=results/hisat2
GTF=raw_files/annotations/turkey_genome_trxpts.gtf
assembled=results/abundances
samples="I_4hrsS1 I_4hrsS2 I_4hrsS3 I_12hrsS1 I_12hrsS3 I_24hrsS1 I_24hrsS2 I_24hrsS3 \
I_72hrsS1 I_72hrsS2 I_72hrsS3 U_4hrsN1 U_4hrsN2 U_12hrsN1 U_12hrsN2 U_24hrsN1 U_24hrsN2 \
U_72hrsN1 U_72hrsN2"

for sample in $samples
do
    echo "Estimating Abundance for ${sample} ..." 
    stringtie -p 10 -eB -G $GTF -o $assembled/abund_${sample}/abund_${sample}.gtf $mapdir/sorted_${sample}.bam
done


secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Transcript Estimations Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs