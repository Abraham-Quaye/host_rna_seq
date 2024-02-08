################# CONVERT .GFF TO GTF AND MOVE AGAT LOGFILE ###########
rule convert_turkey_gff_to_gtf:
    input:
        script = "scripts/shell_code/make_host_gtf.zsh",
        gff = "raw_files/annotations/turkey_genome.gff"
    output:
        "raw_files/annotations/turkey_genome.gtf",
        "raw_files/annotations/turkey_genome.agat.log"
    shell:
        "{input.script}"

################### EXTRACT SPLICE-SITES ####################
rule extract_host_splice_site:
    input:
        script = "scripts/shell_code/extract_ss.zsh",
        gtf = "raw_files/annotations/turkey_genome.gtf"
    output:
        "raw_files/annotations/turkey_genome.ss"
    shell:
        "{input.script}"

################### EXTRACT EXONS  ####################
rule extract_host_exons:
    input:
        script = "scripts/shell_code/extract_exons.zsh",
        gtf = "raw_files/annotations/turkey_genome.gtf"
    output:
        "raw_files/annotations/turkey_genome.exons"
    shell:
        "{input.script}"

#################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
rule build_host_genome_index:
    input:
        script = "scripts/shell_code/build_genome_index.zsh",
        genome = "raw_files/genome_file/turkey_genome.fa",
        exons = "raw_files/annotations/turkey_genome.exons",
        ss = "raw_files/annotations/turkey_genome.ss"
    output:
        expand("raw_files/index/turkey_tran.{num}.ht2", \
        num = range(1, 9))
    shell:
        "{input.script}"

#################### TRIM READS WITH TRIMGALORE #############
rule trim_and_QC_rnaseq_reads:
    input:
        script = "scripts/shell_code/trim_reads.zsh",
        f1 = expand("raw_files/initial_reads/forwardData/I_{tp}hrsS{rep}_Data1.fq.gz", \
        tp = [4, 24, 72], rep = range(1, 4)),
        f2 = expand("raw_files/initial_reads/forwardData/I_12hrsS{rep}_Data1.fq.gz", \
        rep = [1, 3]),
        f3 = expand("raw_files/initial_reads/forwardData/U_{tp}hrsN{rep}_Data1.fq.gz", \
        tp = [4, 12, 24, 72], rep = [1, 2]),
        # data2
        r1 = expand("raw_files/initial_reads/reverseData/I_{tp}hrsS{rep}_Data2.fq.gz", \
        tp = [4, 24, 72], rep = range(1, 4)),
        r2 = expand("raw_files/initial_reads/reverseData/I_12hrsS{rep}_Data2.fq.gz", \
        rep = [1, 3]),
        r3 = expand("raw_files/initial_reads/reverseData/U_{tp}hrsN{rep}_Data2.fq.gz", \
        tp = [4, 12, 24, 72], rep = [1, 2])
    output:
        expand("trimmedReads/I_{tp}hrsS{rep}_Data1_val_1.fq.gz", \
        tp = [4, 24, 72], rep = range(1, 4)),
        expand("trimmedReads/I_12hrsS{rep}_Data1_val_1.fq.gz", \
        rep = [1, 3]),
        expand("trimmedReads/U_{tp}hrsN{rep}_Data1_val_1.fq.gz", \
        tp = [4, 12, 24, 72], rep = [1, 2]),
        # data2
        expand("trimmedReads/I_{tp}hrsS{rep}_Data2_val_2.fq.gz", \
        tp = [4, 24, 72], rep = range(1, 4)),
        expand("trimmedReads/I_12hrsS{rep}_Data2_val_2.fq.gz", \
        rep = [1, 3]),
        expand("trimmedReads/U_{tp}hrsN{rep}_Data2_val_2.fq.gz", \
        tp = [4, 12, 24, 72], rep = [1, 2])
    shell:
        "{input.script}"




#################### MAP READS TO TURKEY GENOME WITH HISAT2 #############
rule map_reads_convert_to_bam:
    input:
        script = "scripts/shell_code/map_sort_to_bam.zsh",
        seqidx = rules.build_host_genome_index.output,
        i_fordata = expand("trimmedReads/I_{tp}hrsS{rep}_Data1_val_1.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        i_for12 = expand("trimmedReads/I_12hrsS{rep}_Data1_val_1.fq.gz", rep = [1, 3]),
        i_revdata = expand("trimmedReads/I_{tp}hrsS{rep}_Data2_val_2.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        i_rev12 = expand("trimmedReads/I_12hrsS{rep}_Data2_val_2.fq.gz", rep = [1, 3]),
        u_fordata = expand("trimmedReads/U_{tp}hrsN{rep}_Data1_val_1.fq.gz", \
        tp = [4, 12, 24, 12], rep = [1, 2]),
        u_revdata = expand("trimmedReads/U_{tp}hrsN{rep}_Data2_val_2.fq.gz", \
        tp = [4, 12, 24, 12], rep = [1, 2])
    output:
        expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam", \
        trtment_tp = ["I_4", "I_24", "I_72"], trt_rep = ["S1", "S2", "S3"]),
        expand("results/hisat2/sorted_I_12hrs{trt_rep}.bam", trt_rep = ["S1", "S3"]),
        expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam", \
        trtment_tp = ["U_4", "U_12", "U_24", "U_72"], trt_rep = ["N1", "N2"])
    shell:
        "{input.script}"

rule run_pipeline:
    input:
        rules.map_reads_convert_to_bam.output

