#!/usr/bin/env bash

f_fastqs=$(ls raw_files/initial_reads/forwardData)
r_fastqs=$(ls raw_files/initial_reads/reverseData)
path=raw_files/initial_reads

echo "sample - reads" > results/fastqc_countF_Reads.txt

for fastq in $f_fastqs;
do
    sample=$(echo $fastq | awk -F '_' '{print($1 "_" $2)}')
    count=$(gunzip --stdout $path/forwardData/${fastq} | grep "@" | wc -l)
    echo "${sample} -${count}" >> results/fastqc_countF_Reads.txt
done

echo "sample - reads" > results/fastqc_countR_Reads.txt

for fastq in $r_fastqs;
do
    sample=$(echo $fastq | awk -F '_' '{print($1 "_" $2)}')
    count=$(gunzip --stdout $path/reverseData/${fastq} | grep "@" | wc -l)
    echo "${sample} -${count}" >> results/fastqc_countR_Reads.txt
done