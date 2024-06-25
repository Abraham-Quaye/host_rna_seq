##########  DOWNLOAD M gallopavo GENOMIC FILES ####################
rule fetch_Mgallopavo_files:
    input:
        "scripts/shell_code/fetch_Mgal_genomeData.zsh"
    output:
        expand("raw_files/annotations/Mgallopavo_{ext}.gz", \
        ext = ["ncbi.gtf", "GOncbi.gaf"]),
        "raw_files/genome_file/Mgallopavo_genome.fna"
    shell:
        "{input}"

########## MODIFY GTF FILE FOR M gallopavo ###############
rule modify_Mgallopavo_GTF:
    input:
        gtf = rules.fetch_Mgallopavo_files.output,
        r_script = "scripts/r_code/modify_ncbi_gtf.R"
    output:
        "raw_files/annotations/mod_turkey_genome.gtf"
    shell:
        "{input.r_script}"

# ##########  REMOVE GENE ATTRIBUTES (EMPTY TRANSCRIPT_ID TAGS) ###########
rule trim_Mgallopavo_gtf_file:
    input:
        rules.modify_Mgallopavo_GTF.output
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
rule extract_splice_sites:
    input:
        rules.modify_Mgallopavo_GTF.output
    output:
        "raw_files/annotations/turkey_genome.ss"
    shell:
        """
        echo "extracting M gallopavo splice-sites ..."
        hisat2_extract_splice_sites.py {input} > {output} &&
        echo "TURKEY splice-site extraction completed successfully"
        """

# ################### EXTRACT EXONS  ####################
rule extract_exons:
    input:
        rules.modify_Mgallopavo_GTF.output
    output:
        "raw_files/annotations/turkey_genome.exons"
    shell:
        """
        echo "M gallopavo extracting exons ..."
        hisat2_extract_exons.py {input} > {output} &&
        echo "TURKEY exon extration complete"
        """

# #################### BUILD THEV GENOMIC INDEX FOR MAPPING WITH HISAT2 ######
rule build_genomic_index:
    input:
        genome = "raw_files/genome_file/Mgallopavo_genome.fna",
        exons = rules.extract_exons.output,
        ss = rules.extract_splice_sites.output
    output:
        expand("raw_files/index/turkey_tran.{num}.ht2", \
        num = range(1, 9))
    # shell:
    #     """
    #     echo "Building TURKEY genomic index ..."
    #     hisat2-build -p 8 --ss {input.ss} --exon {input.exons} {input.genome} raw_files/index/turkey_tran &&
    #     echo "Index built successfully" || echo "Program Aborted!!!"
    #     """

#################### CLEAN FASTQC READS ####################
rule trimmed_fastQ_files:
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

# #################### QC READS WITH FASTQC ####################
rule FastQC_reads:
    input:
       rules.trimmed_fastQ_files.input
    output:
        expand("qc_reports/fastqc_files/I_{tp}hrsS{rep}_Data1_fastqc.{ext}", \
        tp = [4, 24, 72], rep = range(1, 4), ext = ["html", "zip"]),
        expand("qc_reports/fastqc_files/I_12hrsS{rep}_Data1_fastqc.{ext}", \
        rep = [1, 3], ext = ["html", "zip"]),
        expand("qc_reports/fastqc_files/U_{tp}hrsN{rep}_Data1_fastqc.{ext}", \
        tp = [4, 12, 24, 72], rep = [1, 2], ext = ["html", "zip"]),
        # data2
        expand("qc_reports/fastqc_files/I_{tp}hrsS{rep}_Data2_fastqc.{ext}", \
        tp = [4, 24, 72], rep = range(1, 4), ext = ["html", "zip"]),
        expand("qc_reports/fastqc_files/I_12hrsS{rep}_Data2_fastqc.{ext}", \
        rep = [1, 3], ext = ["html", "zip"]),
        expand("qc_reports/fastqc_files/U_{tp}hrsN{rep}_Data2_fastqc.{ext}", \
        tp = [4, 12, 24, 72], rep = [1, 2], ext = ["html", "zip"])
    shell:
        """
        echo "Running FastQC ..."
        fastqc --memory 5120 -t 19 -o qc_reports/fastqc_files {input.f1} {input.f2} {input.f3}
        fastqc --memory 5120 -t 19 -o qc_reports/fastqc_files {input.r1} {input.r2} {input.r3}
        """

#################### QC READS WITH MULTIQC ####################
rule MultiQC_reads:
    input:
        rules.FastQC_reads.output
    output:
        "qc_reports/multiqc_files/multiqc_report.html",
        directory("qc_reports/multiqc_files/multiqc_data")
    shell:
        """
        multiqc -f -o qc_reports/multiqc_files {input}
        """

#################### COUNT TOTAL READS IN FASTQC FILES (TRIMMED READS) #######
rule count_trimmed_reads:
    input:
        fastqc = rules.trimmed_fastQ_files.input,
        script = "scripts/shell_code/fastqc_readCounts.sh"
    output:
        expand("results/fastqc_count{end}_Reads.txt", end = ["F", "R"])
    shell:
        "{input.script}"

# #################### MAP READS TO TURKEY GENOME WITH HISAT2 #############
rule map_reads_convert_to_bam:
    input:
        script = "scripts/shell_code/map_sort_to_bam.zsh",
        seqidx = rules.build_genomic_index.output,
        i_fordata = expand("raw_files/initial_reads/forwardData/I_{tp}hrsS{rep}_Data1.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        i_for12 = expand("raw_files/initial_reads/forwardData/I_12hrsS{rep}_Data1.fq.gz", \
        rep = [1, 3]),
        i_revdata = expand("raw_files/initial_reads/reverseData/I_{tp}hrsS{rep}_Data2.fq.gz", \
        tp = [72, 24, 4], rep = [1, 2, 3]),
        i_rev12 = expand("raw_files/initial_reads/reverseData/I_12hrsS{rep}_Data2.fq.gz", \
        rep = [1, 3]),
        u_fordata = expand("raw_files/initial_reads/forwardData/U_{tp}hrsN{rep}_Data1.fq.gz", \
        tp = [4, 12, 24, 12], rep = [1, 2]),
        u_revdata = expand("raw_files/initial_reads/reverseData/U_{tp}hrsN{rep}_Data2.fq.gz", \
        tp = [4, 12, 24, 12], rep = [1, 2])
    output:
        expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam", \
        trtment_tp = ["I_4", "I_24", "I_72"], trt_rep = ["S1", "S2", "S3"]),
        expand("results/hisat2/sorted_I_12hrs{trt_rep}.bam", trt_rep = ["S1", "S3"]),
        expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam", \
        trtment_tp = ["U_4", "U_12", "U_24", "U_72"], trt_rep = ["N1", "N2"])
    shell:
        "{input.script}"

#################### INDEX SORTED BAM FILES WITH SAMTOOLS #############
rule index_sorted_bamFiles:
    input:
        script = "scripts/shell_code/index_sorted_bams.zsh",
        bams = rules.map_reads_convert_to_bam.output
    output:
        expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam.bai", \
        trtment_tp = ["I_4", "I_24", "I_72"], trt_rep = ["S1", "S2", "S3"]),
        expand("results/hisat2/sorted_I_12hrs{trt_rep}.bam.bai", trt_rep = ["S1", "S3"]),
        expand("results/hisat2/sorted_{trtment_tp}hrs{trt_rep}.bam.bai", \
        trtment_tp = ["U_4", "U_12", "U_24", "U_72"], trt_rep = ["N1", "N2"])
    shell:
        "{input.script}"

#################### MAPPING STATISTICS FOR BAM FILES #######
rule extract_mapping_statistics:
    input:
        fastqc = rules.map_reads_convert_to_bam.output,
        script = "scripts/shell_code/mapping_stats.sh",
        py_script = "scripts/python/calc_GC_content.py"
    output:
        expand("results/count{feat}_Mapped_reads.txt", feat = ["All", "Unq"]),
        expand("results/countQ{q}_reads.txt", q = [20, 30]),
        "results/gc_content_results.csv"
    shell:
        """
        {input.script}
        {input.py_script}
        """

#################### ESTIMATE TRANSCRIPT ABUNDANCES WITH STRINGTIE #############
rule estimate_trancript_abundances:
    input:
        bams = rules.map_reads_convert_to_bam.output,
        merged_gtf = rules.trim_Mgallopavo_gtf_file.output,
        script = "scripts/shell_code/est_trxpt_abund.zsh"
    output:
        expand("results/abundances/abund_{sample}/abund_{sample}.gtf", \
        sample = ["I_4hrsS1", "I_4hrsS2", "I_4hrsS3", "I_12hrsS1", "I_12hrsS3", \
        "I_24hrsS1", "I_24hrsS2", "I_24hrsS3", "I_72hrsS1", "I_72hrsS2", \
        "I_72hrsS3", "U_4hrsN1", "U_4hrsN2", "U_12hrsN1", "U_12hrsN2", \
        "U_24hrsN1", "U_24hrsN2", "U_72hrsN1", "U_72hrsN2"])
    shell:
        "{input.script}"

#################### MAKE COUNT MATRICES WITH STRINGTIE PYTHON SCRIPT #############
rule generate_count_matrices:
    input:
        counts = rules.estimate_trancript_abundances.output,
        prog = "scripts/python/prepDE.py3",
        script = "scripts/shell_code/make_count_matrix.zsh"
    output:
        expand("results/abundances/count_matrix/{feat}_count_matrix.csv", \
        feat = ["genes", "trxpts"])
    shell:
        "{input.script}"

####### EXTRACT AND SAVE GENE IDS FOR GO AND KEGG ANALYSES ###########################
rule DESeq2_DEG_analysis:
    input:
        cnt_matrix = "results/abundances/count_matrix/genes_count_matrix.csv",
        r_script = "scripts/r_code/deseq_analysis.R"
    output:
        sigs = expand("results/r/tables/signif_{tp}hrsDEGs.csv", tp = [4, 12, 24, 72]),
        total = expand("results/r/tables/total_{tp}hrsDEGs.csv", tp = [4, 12, 24, 72]),
        fig = expand("results/r/figures/{type}_{tp}hrs.png",
        type = ["pca", "volcano", "distPlot"], tp = [12, 24]),
        corr_fig = "results/r/figures/sample_corr_figure.png"
    shell:
        """
        {input.r_script}
        """


####### EXTRACT AND SAVE GENE IDS FOR GO AND KEGG ANALYSES ###########################
rule save_DESeq2_results:
    input:
        deg_files = rules.DESeq2_DEG_analysis.output.sigs,
        r_script1 = "scripts/r_code/my_degAnalyses.R",
        r_script2 = "scripts/r_code/extract_myDEG_geneIDs.R",
        r_script3 = "scripts/r_code/save_deg_tabs.R"
    output:
        expand("results/r/tables/{names}.txt", \
        names = ["down_degIDs_4hrs", "down_degIDs_12hrs", "down_degIDs_24hrs", \
        "down_degIDs_72hrs", "up_degIDs_4hrs", "up_degIDs_12hrs", "up_degIDs_24hrs", \
        "all_down_degs", "all_up_degs"])
    shell:
        "{input.r_script3}"

# ####### PLOT COMPOSITE FIGURE FOR UP AND DOWN REGULATED DEGs ###########################
rule plot_DEG_figures:
    input:
        deg_files = rules.DESeq2_DEG_analysis.output.sigs,
        r_script1 = "scripts/r_code/plot_myDEGS.R",
        r_script2 = "scripts/r_code/plot_myHeatmap.R",
        r_script3 = "scripts/r_code/my_degAnalyses.R",
        r_script4 = "scripts/r_code/extract_myDEG_geneIDs.R",
        r_script5 = "scripts/r_code/my_vennDiagram.R"
    output:
        "results/r/figures/deg_patch_fig.png"
    shell:
        """
        {input.r_script2}
        rm Rplots.pdf
        """

####### GO TERM AND PATHWAY ENRICHMENT ANALYSIS ###########################
rule plot_enrichment:
    input:
        r_script1 = "scripts/r_code/extract_myDEG_geneIDs.R",
        r_script2 = "scripts/r_code/my_degAnalyses.R",
        degfiles = rules.DESeq2_DEG_analysis.output.sigs,
        main_script = "scripts/r_code/my_enrichment_analyses.R"
    output:
        go_figs = expand("results/r/figures/go_enrich_{tp}{reg}{GO}.png", \
        tp = [12, 24], reg = ["up", "down"], GO = ["BP", "CC", "MF"]),
        patch_fig = "results/r/figures/patch_GO_enrich.png",
        go_res_tables = expand("results/r/tables/t{tp}_GO_results.tsv", tp = [12, 24])
    shell:
        "{input.main_script}"
    
####### WRITE MANUSCRIPT FOR PUBLICATION ###########################
rule write_manuscript:
    input:
        rmd = "infected_host_trxptome.Rmd",
        ref_style = "asm.csl",
        refs = "trxptome_refs.bib",
        map_stats_script = "scripts/r_code/reads_mapping_stats.R",
        map_stats_data = rules.extract_mapping_statistics.output,
        trimmed_rds = rules.count_trimmed_reads.output,
        fig2 = rules.plot_enrichment.output.patch_fig,
        fig3 = rules.plot_DEG_figures.output,
        go_res = rules.plot_enrichment.output.go_res_tables
    output:
        "infected_host_trxptome.pdf",
        "infected_host_trxptome.tex",
        "infected_host_trxptome.docx"
    shell:
        """
        R -e "library(rmarkdown);render('{input.rmd}', output_format = 'all')";
        """
rule run_pipeline:
    input:
        rules.MultiQC_reads.output,
        rules.save_DESeq2_results.output,
        rules.write_manuscript.output,
        rules.index_sorted_bamFiles.output