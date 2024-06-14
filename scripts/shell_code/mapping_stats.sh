#!/usr/bin/env bash

mapdir=results/hisat2
bams=$(ls results/hisat2/*.bam | awk -F '/' '{print $3}')


# Count all mapped reads
echo "bamfile,mapped_reads" > results/countAll_Mapped_reads.txt
for bam in $bams;
do
    bamfile=${mapdir}/${bam}
    count=$(samtools view -@ 8 -c -F 4 $bamfile) 
    echo "${bam},${count}" >> results/countAll_Mapped_reads.txt
done


# Count only uniquely mapped reads
echo "bamfile,unq_mapped_reads" > results/countUnq_Mapped_reads.txt
for bam in $bams;
do
    bamfile=${mapdir}/${bam}
    count=$(samtools view -@ 8 -c -F 260 $bamfile) 
    echo "${bam},${count}" >> results/countUnq_Mapped_reads.txt
done


# Count reads with mapping quality of at least 30
echo "bamfile,q30" > results/countQ30_reads.txt
for bam in $bams;
do
    bamfile=${mapdir}/${bam}
    count=$(samtools view -@ 8 -c -q 30 $bamfile) 
    echo "${bam},${count}" >> results/countQ30_reads.txt
done

# Count reads with mapping quality of at least 30
echo "bamfile,q20" > results/countQ20_reads.txt
for bam in $bams;
do
    bamfile=${mapdir}/${bam}
    count=$(samtools view -@ 8 -c -q 20 $bamfile) 
    echo "${bam},${count}" >> results/countQ20_reads.txt
done