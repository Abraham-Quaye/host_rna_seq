#!/usr/bin/env zsh

seqidx=raw_files/index/turkey_tran
mapdir=results/hisat2

# paired-end reads data1
forData=(I_4hrsS1_Data1_val_1.fq.gz I_4hrsS2_Data1_val_1.fq.gz I_4hrsS3_Data1_val_1.fq.gz \
I_12hrsS1_Data1_val_1.fq.gz I_12hrsS3_Data1_val_1.fq.gz \
I_24hrsS1_Data1_val_1.fq.gz I_24hrsS2_Data1_val_1.fq.gz I_24hrsS3_Data1_val_1.fq.gz \
I_72hrsS1_Data1_val_1.fq.gz I_72hrsS2_Data1_val_1.fq.gz I_72hrsS3_Data1_val_1.fq.gz \
U_4hrsN1_Data1_val_1.fq.gz U_4hrsN2_Data1_val_1.fq.gz \
U_12hrsN1_Data1_val_1.fq.gz U_12hrsN2_Data1_val_1.fq.gz \
U_24hrsN1_Data1_val_1.fq.gz U_24hrsN2_Data1_val_1.fq.gz \
U_72hrsN1_Data1_val_1.fq.gz U_72hrsN2_Data1_val_1.fq.gz)

# paired-end reads data2
revData=(I_4hrsS1_Data2_val_2.fq.gz I_4hrsS2_Data2_val_2.fq.gz I_4hrsS3_Data2_val_2.fq.gz \
I_12hrsS1_Data2_val_2.fq.gz I_12hrsS3_Data2_val_2.fq.gz \
I_24hrsS1_Data2_val_2.fq.gz I_24hrsS2_Data2_val_2.fq.gz I_24hrsS3_Data2_val_2.fq.gz \
I_72hrsS1_Data2_val_2.fq.gz I_72hrsS2_Data2_val_2.fq.gz I_72hrsS3_Data2_val_2.fq.gz \
U_4hrsN1_Data2_val_2.fq.gz U_4hrsN2_Data2_val_2.fq.gz \
U_12hrsN1_Data2_val_2.fq.gz U_12hrsN2_Data2_val_2.fq.gz \
U_24hrsN1_Data2_val_2.fq.gz U_24hrsN2_Data2_val_2.fq.gz \
U_72hrsN1_Data2_val_2.fq.gz U_72hrsN2_Data2_val_2.fq.gz)

# arrays for labelling output files
trt_tps=(I_4 I_4 I_4 I_12 I_12 I_24 I_24 I_24 I_72 I_72 I_72 U_4 U_4 U_12 U_12 U_24 U_24 U_72 U_72)
lab_reps=(S1 S2 S3 S1 S3 S1 S2 S3 S1 S2 S3 N1 N2 N1 N2 N1 N2 N1 N2)

for (( s = 1; s <=19; s++ ))
do
    name=${trt_tps[$s]}
    lb_rep=${lab_reps[$s]}
    d1=${forData[$s]}
    d2=${revData[$s]}
    echo "Mapping Reads to ${name}hrs${lb_rep} Now..."
    
    echo "$d1 goes with $d2 -- named ${name}${lb_rep}"
    # hisat2 -p 10 --dta -x $seqidx -1 trimmedReads/$d1 -2 trimmedReads/$d2 \
    # | samtools sort -@ 10 -o $mapdir/sorted_${name}hrs${lb_rep}.bam

        secs=$SECONDS
        hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
        echo "Completed Mapping to ${name}${lb_rep} in: ${hrs}:${mins}:${secs}"
done

secs=$SECONDS
hrs=$(( secs/3600 )); mins=$(( (secs-hrs*3600)/60 )); secs=$(( secs-hrs*3600-mins*60 ))
printf 'Read Mapping Script Completed in: %02d:%02d:%02d\n' $hrs $mins $secs