#!/bin/bash

forReads=raw_files/initial_reads/forwardData
revReads=raw_files/initial_reads/reverseData

function shave_reads(){
    echo "Running TrimGalore!..."
    read1=$forReads/$1
    read2=$revReads/$2
    trim_galore --phred33 -q 20 --fastqc --gzip --cores 6 --paired $read1 $read2 -o trimmedReads/
}

echo "Trimming 4hr Samples"
shave_reads I_4hrsS1_Data1.fq.gz I_4hrsS1_Data2.fq.gz
shave_reads I_4hrsS2_Data1.fq.gz I_4hrsS2_Data2.fq.gz 
shave_reads I_4hrsS3_Data1.fq.gz I_4hrsS3_Data2.fq.gz 
# uninfected
shave_reads U_4hrsN1_Data1.fq.gz U_4hrsN1_Data2.fq.gz
shave_reads U_4hrsN2_Data1.fq.gz U_4hrsN2_Data2.fq.gz

echo "Trimming 12hr Samples"
shave_reads I_12hrsS1_Data1.fq.gz I_12hrsS1_Data2.fq.gz
shave_reads I_12hrsS3_Data1.fq.gz I_12hrsS3_Data2.fq.gz 
# uninfected
shave_reads U_12hrsN1_Data1.fq.gz U_12hrsN1_Data2.fq.gz
shave_reads U_12hrsN2_Data1.fq.gz U_12hrsN2_Data2.fq.gz

echo "Trimming 24hr Samples"
shave_reads I_24hrsS1_Data1.fq.gz I_24hrsS1_Data2.fq.gz
shave_reads I_24hrsS2_Data1.fq.gz I_24hrsS2_Data2.fq.gz 
shave_reads I_24hrsS3_Data1.fq.gz I_24hrsS3_Data2.fq.gz 
# uninfected
shave_reads U_24hrsN1_Data1.fq.gz U_24hrsN1_Data2.fq.gz
shave_reads U_24hrsN2_Data1.fq.gz U_24hrsN2_Data2.fq.gz

echo "Trimming 72hr Samples"
shave_reads I_72hrsS1_Data1.fq.gz I_72hrsS1_Data2.fq.gz
shave_reads I_72hrsS2_Data1.fq.gz I_72hrsS2_Data2.fq.gz 
shave_reads I_72hrsS3_Data1.fq.gz I_72hrsS3_Data2.fq.gz 
# uninfected
shave_reads U_72hrsN1_Data1.fq.gz U_72hrsN1_Data2.fq.gz
shave_reads U_72hrsN2_Data1.fq.gz U_72hrsN2_Data2.fq.gz

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Script completed in: %02d:%02d:%02d\n' $hrs $mins $secs
