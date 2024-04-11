##########  DOWNLOAD M gallopavo GENOMIC FILES ####################
rule fetch_Mgal_files:
    input:
        "scripts/shell_code/fetch_Mgal_genomeData.zsh"
    output:
        expand("raw_files/annotations/Mgallopavo_{ext}.gz", \
        ext = ["ncbi.gtf", "GOncbi.gaf"]),
        "raw_files/genome_file/Mgallopavo_genome.fna"
    shell:
        "{input}"

##########  MAKE LOCAL ORGANISM DATABASE FOR M gallopavo ###############
rule make_Mgallopavo_OrgDB:
    input:
        data = "raw_analysis/count_matrix/all_genes_fpkm_exp.xlsx",
        r_script = "scripts/r_code/install_Mga_annotationDB.R"
    output:
       directory("org.Mgallopavo.eg.db")
    shell:
        "{input.r_script}"

########## MODIFY GTF FILE FOR M gallopavo ###############
rule modify_turkey_GTF:
    input:
        gtf = rules.fetch_Mgal_files.output,
        r_script = "scripts/r_code/modify_ncbi_gtf.R"
    output:
        "raw_files/annotations/mod_turkey_genome.gtf"
    shell:
        "{input.r_script}"

# ##########  REMOVE GENE ATTRIBUTES (EMPTY TRANSCRIPT_ID TAGS) ###########
rule trim_turkey_gtf_file:
    input:
        rules.modify_turkey_GTF.output
    output:
        "raw_files/annotations/turkey_genome_trxpts.gtf"
    shell:
        """
        echo "Creating Trimmed GTF file"
        awk '$3 != "gene" ' {input} > {output} &&
        echo "GTF with no empty tags created successfully!!!" \
        || echo "Error in creating trimmed GTF file"
        """

# ################### EXTRACT SPLICE-SITES ####################
rule extract_host_splice_site:
    input:
        rules.modify_turkey_GTF.output
    output:
        "raw_files/annotations/turkey_genome.ss"
    shell:
        """
        echo "extracting M gallopavo splice-sites ..."
        hisat2_extract_splice_sites.py {input} > {output} &&
        echo "TURKEY splice-site extraction completed successfully"
        """

# ################### EXTRACT EXONS  ####################
rule extract_host_exons:
    input:
        rules.modify_turkey_GTF.output
    output:
        "raw_files/annotations/turkey_genome.exons"
    shell:
        """
        echo "M gallopavo extracting exons ..."
        hisat2_extract_exons.py {input} > {output} &&
        echo "TURKEY exon extration complete"
        """

# #################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
rule build_host_genome_index:
    input:
        genome = "raw_files/genome_file/Mgallopavo_genome.fna",
        exons = rules.extract_host_exons.output,
        ss = rules.extract_host_splice_site.output
    output:
        expand("raw_files/index/turkey_tran.{num}.ht2", \
        num = range(1, 9))
    shell:
        """
        echo "Building TURKEY genomic index ..."
        hisat2-build -p 8 --ss {input.ss} --exon {input.exons} {input.genome} raw_files/index/turkey_tran &&
        echo "Index built successfully" || echo "Program Aborted!!!"
        """

# #################### TRIM READS WITH TRIMGALORE #############
rule FastQC_reads:
    input:
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
        expand("fastqc_files/I_{tp}hrsS{rep}_Data1_fastqc.{ext}", \
        tp = [4, 24, 72], rep = range(1, 4), ext = ["html", "zip"]),
        expand("fastqc_files/I_12hrsS{rep}_Data1_fastqc.{ext}", \
        rep = [1, 3], ext = ["html", "zip"]),
        expand("fastqc_files/U_{tp}hrsN{rep}_Data1_fastqc.{ext}", \
        tp = [4, 12, 24, 72], rep = [1, 2], ext = ["html", "zip"]),
        # data2
        expand("fastqc_files/I_{tp}hrsS{rep}_Data2_fastqc.{ext}", \
        tp = [4, 24, 72], rep = range(1, 4), ext = ["html", "zip"]),
        expand("fastqc_files/I_12hrsS{rep}_Data2_fastqc.{ext}", \
        rep = [1, 3], ext = ["html", "zip"]),
        expand("fastqc_files/U_{tp}hrsN{rep}_Data2_fastqc.{ext}", \
        tp = [4, 12, 24, 72], rep = [1, 2], ext = ["html", "zip"])
    shell:
        """
        echo "Running FastQC ..."
        fastqc --memory 5120 -t 19 -o fastqc_files {input.f1} {input.f2} {input.f3}
        fastqc --memory 5120 -t 19 -o fastqc_files {input.r1} {input.r2} {input.r3}
        """

# #################### MAP READS TO TURKEY GENOME WITH HISAT2 #############
# rule map_reads_convert_to_bam:
#     input:
#         script = "scripts/shell_code/map_sort_to_bam.zsh",
#         seqidx = rules.build_host_genome_index.output,
#         i_fordata = expand("trimmedReads/I_{tp}hrsS{rep}_Data1_val_1.fq.gz", \
#         tp = [72, 24, 4], rep = [1, 2, 3]),
#         i_for12 = expand("trimmedReads/I_12hrsS{rep}_Data1_val_1.fq.gz", rep = [1, 3]),
#         i_revdata = expand("trimmedReads/I_{tp}hrsS{rep}_Data2_val_2.fq.gz", \
#         tp = [72, 24, 4], rep = [1, 2, 3]),
#         i_rev12 = expand("trimmedReads/I_12hrsS{rep}_Data2_val_2.fq.gz", rep = [1, 3]),
#         u_fordata = expand("trimmedReads/U_{tp}hrsN{rep}_Data1_val_1.fq.gz", \
#         tp = [4, 12, 24, 12], rep = [1, 2]),
#         u_revdata = expand("trimmedReads/U_{tp}hrsN{rep}_Data2_val_2.fq.gz", \
#         tp = [4, 12, 24, 12], rep = [1, 2])
#     output:
#         expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam", \
#         trtment_tp = ["I_4", "I_24", "I_72"], trt_rep = ["S1", "S2", "S3"]),
#         expand("results/hisat2/sorted_I_12hrs{trt_rep}.bam", trt_rep = ["S1", "S3"]),
#         expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam", \
#         trtment_tp = ["U_4", "U_12", "U_24", "U_72"], trt_rep = ["N1", "N2"])
#     shell:
#         "{input.script}"

# #################### INDEX SORTED BAM FILES WITH SAMTOOLS #############
# rule index_sorted_bamFiles:
#     input:
#         script = "scripts/shell_code/index_sorted_bams.zsh",
#         bams = rules.map_reads_convert_to_bam.output
#     output:
#         expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam.bai", \
#         trtment_tp = ["I_4", "I_24", "I_72"], trt_rep = ["S1", "S2", "S3"]),
#         expand("results/hisat2/sorted_I_12hrs{trt_rep}.bam.bai", trt_rep = ["S1", "S3"]),
#         expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam.bai", \
#         trtment_tp = ["U_4", "U_12", "U_24", "U_72"], trt_rep = ["N1", "N2"])
#     shell:
#         "{input.script}"

# #################### ASSEMBLE TRANSCRIPTS WITH STRINGTIE #############
# rule assemble_transcripts:
#     input:
#         bams = rules.map_reads_convert_to_bam.output,
#         gtf = rules.trim_turkey_gtf_file.output,
#         script = "scripts/shell_code/assemble_transcripts.zsh"
#     output:
#         expand("results/stringtie/I_{tp}hrsS{rep}.gtf", tp = [4, 24, 72], rep = range(1, 4)),
#         expand("results/stringtie/I_12hrsS{rep}.gtf", rep = [1, 3]),
#         expand("results/stringtie/U_{tp}hrsN{rep}.gtf", tp = [4, 12, 24, 72], rep = [1, 2])
#     shell:
#         "{input.script}"

# #################### LIST ALL STRINGTIE GTF FILES TO BE MERGED MERGE #############
# rule transcript_merge_list:
#     input:
#         gtfs_sngles = rules.assemble_transcripts.output
#     output:
#         "results/stringtie/trxpt_merge_list.txt"
#     shell:
#         """
#         ls results/stringtie/*.gtf > {output}
#         """

# #################### MERGE ALL TRANSCRIPTS WITH STRINGTIE #############
# rule merge_assembled_transcripts:
#     input:
#         gtfs_sngles = rules.assemble_transcripts.output,
#         gtf_list = rules.transcript_merge_list.output,
#         gtf_main = rules.trim_turkey_gtf_file.output
#     output:
#         "results/stringtie/turkey_merged_all_tps.gtf"
#     shell:
#         """
#         echo "Merging all transcripts from all timepoints now...";
#         stringtie --merge -p 8 -G {input.gtf_main} -l rp19 -o {output} {input.gtf_list};
#         echo "Merging Completed in $SECONDS secs"
#         """

# #################### COMPARE MERGED TRANSCRIPTS WITH REFERENCE TRANSCRIPTS #############
# rule compare_merged_trxpts_toReference:
#     input:
#         gtf_merged = rules.merge_assembled_transcripts.output,
#         gtf_ref = rules.trim_turkey_gtf_file.output
#     output:
#         expand("results/gffcompare/turkey_merged.{type}", \
#         type = ["stats", "loci", "tracking", "annotated.gtf", \
#         "turkey_merged_all_tps.gtf.refmap", "turkey_merged_all_tps.gtf.tmap"])
#     shell:
#         """
#         echo "Comparing merged transcripts with reference...";
#         # -G flag = tells gffcompare to compare all transcripts in the input transcripts.gtf file

#         gffcompare -G -r {input.gtf_ref} -o turkey_merged {input.gtf_merged}
#         mv turkey_merged* results/gffcompare;
#         mv results/stringtie/turkey_merged.turkey_merged_all_tps* results/gffcompare;
#         echo "gffcompare completed in $SECONDS secs"
#         """

# #################### ESTIMATE TRANSCRIPT ABUNDANCES WITH STRINGTIE #############
# rule estimate_trancript_abundances:
#     input:
#         bams = rules.map_reads_convert_to_bam.output,
#         merged_gtf = rules.merge_assembled_transcripts.output,
#         script = "scripts/shell_code/est_trxpt_abund.zsh"
#     output:
#         expand("results/abundances/abund_{sample}/abund_{sample}.gtf", \
#         sample = ["I_4hrsS1", "I_4hrsS2", "I_4hrsS3", "I_12hrsS1", "I_12hrsS3", \
#         "I_24hrsS1", "I_24hrsS2", "I_24hrsS3", "I_72hrsS1", "I_72hrsS2", \
#         "I_72hrsS3", "U_4hrsN1", "U_4hrsN2", "U_12hrsN1", "U_12hrsN2", \
#         "U_24hrsN1", "U_24hrsN2", "U_72hrsN1", "U_72hrsN2"])
#     shell:
#         "{input.script}"

# #################### MAKE COUNT MATRICES WITH STRINGTIE PYTHON SCRIPT #############
# rule generate_count_matrices:
#     input:
#         counts = rules.estimate_trancript_abundances.output,
#         prog = "scripts/python/prepDE.py3",
#         script = "scripts/shell_code/make_count_matrix.zsh"
#     output:
#         expand("results/abundances/count_matrix/{feat}_count_matrix.csv", \
#         feat = ["genes", "trxpts"])
#     shell:
#         "{input.script}"

# ####### PLOTS FOR DIFFERENTIAL GENE EXPRESSION ################################
# rule plot_DEG_and_Heatmaps:
#     input:
#         counts = rules.generate_count_matrices.output, # tentative connection
#         deg_files = expand("raw_analysis/diff_exp/diff_gene_exp_{tp}hrs.xlsx", \
#         tp = [4, 12, 24, 72]),
#         r_script1 = "scripts/r_code/plot_degs.R",
#         r_script2 = "scripts/r_code/plot_heatmaps.R"
#     output:
#         "results/r/figures/deg_bar_plt.png",
#         expand("results/r/figures/deg_heatmap_{tpp}hpi.png", tpp = [4, 12, 24])
#     shell:
#         """
#         {input.r_script1}
#         {input.r_script2}
#         rm Rplots.pdf
#         """

# ####### EXTRACT AND SAVE GENE IDS FOR GO AND KEGG ANALYSES ###########################
# rule extract_geneIDs_GO_KEGG:
#     input:
#         counts = rules.generate_count_matrices.output,
#         r_script1 = "scripts/r_code/deg_analysis.R",
#         deg_files = expand("raw_analysis/diff_exp/diff_gene_exp_{tp}hrs.xlsx", \
#         tp = [4, 12, 24, 72]),
#         r_script2 = "scripts/r_code/extract_deg_geneIDs.R",
#         r_script3 = "scripts/r_code/save_deg_tabs.R"
#     output:
#         expand("results/r/tables/{names}.txt", \
#         names = ["all_bg_geneIDs", "down_degIDs_4hrs", "down_degIDs_12hrs", \
#         "up_degIDs_12hrs", "down_degIDs_24hrs", "up_degIDs_24hrs", "all_down_degs", \
#         "all_up_degs"])
#     shell:
#         "{input.r_script3}"

rule run_pipeline:
    input:
        rules.make_Mgallopavo_OrgDB.output,
        rules.trim_turkey_gtf_file.output,
        rules.build_host_genome_index.output,
        rules.FastQC_reads.output
#         rules.compare_merged_trxpts_toReference.output,
#         rules.index_sorted_bamFiles.output,
#         rules.plot_DEG_and_Heatmaps.output,
#         rules.extract_geneIDs_GO_KEGG.output

