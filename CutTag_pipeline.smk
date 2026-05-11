# Snakefile

# By Goboru, November 2025
# Pipeline to run the complete pipeline to analyze ATAC data or do individual steps 


# Use: conda activate ATACpipeline; snakemake -s ATAC_pipeline.smk --cores 8

# Configurable paths
configfile: "config.yaml"  # optional

# CONFIGURATION
dir_raw = config.get("raw_dir", "raw_data")
dir_out = config.get("out_dir", "results")

# Gather FASTQ files
SAMPLES, = glob_wildcards(f"{dir_raw}/{{sample}}.fastq.gz")
UNIQ_SAMPLES, = glob_wildcards(f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz")

wildcard_constraints:
    uniq_sample = "((?!_r1).)*"

# This rule sets when the pipeline is finished
rule all:
    input: 
        f"{dir_out}/qc_untrimmed/multiqc_report.html",
        f"{dir_out}/qc_trimmed/multiqc_report.html",
        expand(f"{dir_out}/no_blacklist/{{uniq_sample}}.noMT.noBlacklist.bam", uniq_sample=UNIQ_SAMPLES),
        expand(f"{dir_out}/no_duplicates/{{uniq_sample}}.final.dedup.bam", uniq_sample=UNIQ_SAMPLES),
        expand(f"{dir_out}/idx_report/{{uniq_sample}}.idxstats.txt", uniq_sample=UNIQ_SAMPLES),
        expand(f"{dir_out}/final_bam_report/{{uniq_sample}}.idxstats.txt", uniq_sample=UNIQ_SAMPLES),
        expand(f"{dir_out}/final_bam_report/{{uniq_sample}}.flagstat.txt", uniq_sample=UNIQ_SAMPLES),
        expand(f"{dir_out}/peaks/{{uniq_sample}}_peaks.narrowPeak", uniq_sample=UNIQ_SAMPLES),
        expand(f"{dir_out}/frip/{{uniq_sample}}_frip.txt", uniq_sample=UNIQ_SAMPLES), 
        f"{dir_out}/frip/all_samples_frip_mqc.png",
        expand(f"{dir_out}/tss/{{uniq_sample}}_tss_enrichment.png", uniq_sample=UNIQ_SAMPLES),
        f"{dir_out}/plots/all_samples_tss_rich.png"



# Rule 1: FastQC before trimming
rule fastqc_untrimmed:
    input:
        f"{dir_raw}/{{sample}}.fastq.gz"
    output:
        html=f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.html",
        zip=f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.zip"
    log:
        f"{dir_out}/logs/qc_untrimmed/{{sample}}.log"
    shell:
        "fastqc {input} --outdir {dir_out}/qc_untrimmed &> {log}"

# Rule 1.2: multiqc of FastQC's before trimming
rule multiqc_untrimmed:
    input:
        expand(f"{dir_out}/qc_untrimmed/{{sample}}_fastqc.zip", sample=SAMPLES)
    output:
        html=f"{dir_out}/qc_untrimmed/multiqc_report.html"
    shell:
        "multiqc {dir_out}/qc_untrimmed --force -o {dir_out}/qc_untrimmed"




# Rule 2: Trimming
rule trimming:
    input:
        r1 = f"{dir_raw}/{{uniq_sample}}_r1.fastq.gz",
        r2 = f"{dir_raw}/{{uniq_sample}}_r2.fastq.gz"
    params:
        mode = config["trim_mode"]
    output:
        r1 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r1.trimmed.fastq.gz",
        r2 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r2.trimmed.fastq.gz",
        json = f"{dir_out}/temp_trimming/{{uniq_sample}}.fastp.json",
        html = f"{dir_out}/temp_trimming/{{uniq_sample}}.fastp.html"
    log:
        f"{dir_out}/logs/trimming/{{uniq_sample}}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            {params.mode} \
            -j {output.json} -h {output.html} \
            &> {log}
        """


# Rule 3: FastQC after trimming
rule fastqc_trimmed:
    input:
        f"{dir_out}/temp_trimming/{{sample}}.trimmed.fastq.gz"
    output:
        html=f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.html",
        zip=f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.zip"
    log:
        f"{dir_out}/logs/qc_trimmed/{{sample}}.log"
    shell:
        "fastqc {input} --outdir {dir_out}/qc_trimmed &> {log}"

# Rule 3.2: multiqc of FastQC's after trimming
rule multiqc_trimmed:
    input:
        expand(f"{dir_out}/qc_trimmed/{{sample}}.trimmed_fastqc.zip", sample=SAMPLES)
    output:
         html=f"{dir_out}/qc_trimmed/multiqc_report.html"
    shell:
        "multiqc {dir_out}/qc_trimmed --force -o {dir_out}/qc_trimmed"


# Rule 4: Alignment with bowtie2

rule bowtie2_align:
    input:
        r1 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r1.trimmed.fastq.gz",
        r2 = f"{dir_out}/temp_trimming/{{uniq_sample}}_r2.trimmed.fastq.gz"
    output:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    log:
        bowtie = f"{dir_out}/logs/bowtie2/{{uniq_sample}}.log",
        samtools = f"{dir_out}/logs/samtools_sort/{{uniq_sample}}.log"
    params:
        index = config["bowtie2_index"]
    threads: config["bowtie2_threads"]
    shell:
        """
        (bowtie2 --very-sensitive -I 25 -X 700 \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} \
            2> {log.bowtie}) \
        | samtools sort \
            -@ 2 \
            -o {output.bam} \
            - \
            2> {log.samtools}
        """
        




# section 5: Post-alingments QCs

#Rule 5.1: idxstats reports
rule idxstats:
    input:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    output:
        report = f"{dir_out}/idx_report/{{uniq_sample}}.idxstats.txt"
    log:
        f"{dir_out}/logs/idx_report/{{uniq_sample}}.idxstats.log"
    threads: 1
    shell:
        """
        # Make sure the BAM is indexed
        samtools index {input.bam} 2>> {log}

        # Generate idxstats report
        samtools idxstats {input.bam} > {output.report} 2>> {log}
        """

# Rule 5.1.1: Plot Mithocondrial Percentage
rule plot_mit_perc:
    input:
        # We still need this expand so Snakemake knows to finish all samples first
        mit_files = expand(f"{dir_out}/idx_report/{{uniq_sample}}.idxstats.txt", uniq_sample=UNIQ_SAMPLES)
    output:
        plot = f"{dir_out}/plots/all_samples_mit_perc.png"
    log:
        f"{dir_out}/logs/plots/mit_perc_plot.log"
    params:
        mit_perc = f"{dir_out}/final_bam_report"
    shell:
        """
        # We pass the directory path (params.frip_dir) instead of the list
        Rscript mit_perc_plot.R {output.plot} {params.mit_perc} &> {log}
        """
# Rule 5.2: remove mitochondrial reads
rule remove_mito:
    input:
        bam = f"{dir_out}/aligned/{{uniq_sample}}_align.bam"
    output:
        bam = f"{dir_out}/mito/{{uniq_sample}}.noMT.bam"
    log:
        f"{dir_out}/logs/mito/{{uniq_sample}}.remove_mito.log"
    threads: 4
    shell:
        r"""
        samtools view -h {input.bam} \
            | awk '$3 != "chrM" && $3 != "MT" || $1 ~ /^@/' \
            | samtools view -b -o {output.bam} \
            &> {log}
        samtools index {output.bam}
        """

# Rule 5.3: Remove ENCODE blacklist regions
rule remove_blacklist:
    input:
        bam = f"{dir_out}/mito/{{uniq_sample}}.noMT.bam"
    output:
        bam = f"{dir_out}/no_blacklist/{{uniq_sample}}.noMT.noBlacklist.bam"
    log:
        f"{dir_out}/logs/no_blacklist/{{uniq_sample}}.remove_blacklist.log"
    params:
        blacklist = config["blacklist_bed"]
    threads: 4
    shell:
        r"""
        bedtools intersect \
            -v \
            -abam {input.bam} \
            -b {params.blacklist} \
            > {output.bam} 2> {log}

        samtools index {output.bam}
        """
# Remove non canonical reads? => Before peak calling


# # Rule 5.4: Mark duplicates, remove duplicates
# rule mark_duplicates:
#     input:
#         bam = f"{dir_out}/no_blacklist/{{uniq_sample}}.noMT.noBlacklist.bam"
#     output:
#         bam  = f"{dir_out}/no_duplicates/{{uniq_sample}}.final.dedup.bam",
#         metrics = f"{dir_out}/no_duplicates/{{uniq_sample}}.final.dedup.metrics.txt"
#     log:
#         f"{dir_out}/logs/no_duplicates/{{uniq_sample}}.dedup.log"
#     threads: 4
#     shell:
#         r"""
#         (picard MarkDuplicates \
#             I={input.bam} \
#             O={output.bam} \
#             M={output.metrics} \
#             REMOVE_DUPLICATES=true \
#             CREATE_INDEX=true  \
#             2> {log})
#         """



#Changing to samtools markdup. requires an extra step to know the duplicate rate

# Rule 5.4: Mark and Remove duplicates using Samtools
rule samtools_dedup:
    input:
        bam = f"{dir_out}/no_blacklist/{{uniq_sample}}.noMT.noBlacklist.bam"
    output:
        bam = f"{dir_out}/no_duplicates/{{uniq_sample}}.final.dedup.bam"
    log:
        f"{dir_out}/logs/no_duplicates/{{uniq_sample}}.sam_dedup.log"
    threads: 8
    shell:
        """
        # 1. Sort by name (required for fixmate)
        samtools sort -n -@ {threads} {input.bam} -o {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.namesort.bam 2> {log}

        # 2. Add mate tags (required for markdup)
        samtools fixmate -m {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.namesort.bam {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.fixmate.bam 2>> {log}

        # 3. Sort by coordinates again
        samtools sort -@ {threads} {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.fixmate.bam -o {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.coordsort.bam 2>> {log}

        # 4. Mark and remove duplicates (-r flag removes them)
        samtools markdup -r -@ {threads} {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.coordsort.bam {output.bam} 2>> {log}

        # 5. Index the final BAM
        samtools index {output.bam} 2>> {log}

        # 6. Cleanup temporary files
        rm {dir_out}/no_duplicates/{{wildcards.uniq_sample}}.tmp.*.bam
        """

# Rule 5.5: Final QC report before peak calling
rule final_bam_qc:
    input:
        bam = f"{dir_out}/no_duplicates/{{uniq_sample}}.final.dedup.bam"
    output:
        idxstats = f"{dir_out}/final_bam_report/{{uniq_sample}}.idxstats.txt",
        flagstat = f"{dir_out}/final_bam_report/{{uniq_sample}}.flagstat.txt"
    log:
        f"{dir_out}/logs/final_bam_report/{{uniq_sample}}.final_bam_qc.log"
    threads: 1
    shell:
        r"""
        # Index BAM (required for idxstats)
        samtools index {input.bam} 2>> {log}

        # idxstats
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}

        # flagstat
        samtools flagstat {input.bam} > {output.flagstat} 2>> {log}
        """

# Rule 6.1: Prepare BAM for MACS2 (Filtering and Sorting)
rule filter_for_macs2:
    input:
        bam = f"{dir_out}/no_duplicates/{{uniq_sample}}.final.dedup.bam"
    output:
        bam = f"{dir_out}/macs2_input/{{uniq_sample}}.filtered.bam"
    log:
        f"{dir_out}/logs/macs2_filter/{{uniq_sample}}.log"
    threads: 4
    shell:
        """
        # Define the primary chromosomes we want to keep
        CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

        # -h: include header
        # -q 30: Quality score >= 30
        # -f 2: Proper pairs only
        # -F 1804: Exclude unmapped, secondary, QC fail, and dups
        # Filter for quality, proper pairs, and specific chromosomes
        samtools view -h -q 30 -f 2 -F 1804 {input.bam} $CHRS | \
        samtools sort -@ {threads} -O bam -o {output.bam} - 
        
        samtools index {output.bam}
        """

# Rule 6.2: Peak calling with MACS2 (Using the filtered BAM)
rule macs2_peak_calling:
    input:
        bam = f"{dir_out}/macs2_input/{{uniq_sample}}.filtered.bam"
    output:
        peaks = f"{dir_out}/peaks/{{uniq_sample}}_peaks.narrowPeak",
        summits = f"{dir_out}/peaks/{{uniq_sample}}_summits.bed",
        xls = f"{dir_out}/peaks/{{uniq_sample}}_peaks.xls"
    log:
        f"{dir_out}/logs/macs2/{{uniq_sample}}.log"
    params:
        out_dir = f"{dir_out}/peaks",
        name = "{uniq_sample}"
    shell:
        """
        macs3 callpeak \
            -t {input.bam} \
            -f BAMPE \
            -g "hs" \
            -n {params.name} \
            --outdir {params.out_dir} \
            -q 0.05 \
            --keep-dup all \
            &> {log}
        """


# Rule 7: Quality controls after peak calling: FRiP 
# Rule 7.1: Calculate FRiP Score
rule calculate_frip:
    input:
        bam = f"{dir_out}/macs2_input/{{uniq_sample}}.filtered.bam",
        peaks = f"{dir_out}/peaks/{{uniq_sample}}_peaks.narrowPeak"
    output:
        frip = f"{dir_out}/frip/{{uniq_sample}}_frip.txt"
    log:
        f"{dir_out}/logs/frip/{{uniq_sample}}.log"
    threads: 4
    shell:
        """
        # 1. Convert narrowPeak to a SAF format (required by featureCounts)
        # SAF format: GeneID  Chr  Start  End  Strand
        awk 'BEGIN{{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"}} \
            {{print $4, $1, $2+1, $3, "."}}' {input.peaks} > {input.peaks}.saf

        # 2. Count reads in peaks using featureCounts
        # -p: paired-end
        # -F SAF: input is in SAF format
        # -a: the annotation (our peaks)
        featureCounts -p -T {threads} -F SAF -a {input.peaks}.saf \
            -o {input.bam}.featureCounts.txt {input.bam} &> {log}

        # 3. Extract the counts and calculate the ratio
        # Total reads is the sum of mapped reads from flagstat or the featureCounts summary
        # We'll use the summary file generated by featureCounts
        READS_IN_PEAKS=$(awk 'NR>2 {{sum+=$7}} END {{print sum}}' {input.bam}.featureCounts.txt)
        TOTAL_READS=$(samtools view -c {input.bam})
        
        # Calculate FRiP
        python3 -c "print(f'{{$READS_IN_PEAKS / $TOTAL_READS:.4f}}')" > {output.frip}

        # 4. Clean up
        rm {input.peaks}.saf {input.bam}.featureCounts.txt {input.bam}.featureCounts.txt.summary
        """

# Rule 7.2: Aggregate FRiP scores and plot
rule plot_all_frip:
    input:
        # We still need this expand so Snakemake knows to finish all samples first
        frip_files = expand(f"{dir_out}/frip/{{uniq_sample}}_frip.txt", uniq_sample=UNIQ_SAMPLES)
    output:
        plot = f"{dir_out}/frip/all_samples_frip_mqc.png"
    log:
        f"{dir_out}/logs/plots/frip_plot.log"
    params:
        frip_dir = f"{dir_out}/frip"
    shell:
        """
        # We pass the directory path (params.frip_dir) instead of the list
        Rscript plot_frip.R {output.plot} {params.frip_dir} &> {log}
        """




# Rule 7.3: Change to bigwig. (Needed for 7.4)
rule bam_to_bw:
    input:
        bam = f"{dir_out}/macs2_input/{{uniq_sample}}.filtered.bam",
        bai = f"{dir_out}/macs2_input/{{uniq_sample}}.filtered.bam.bai"
    output:
        bw = f"{dir_out}/bigwig/{{uniq_sample}}.bw"
    log:
        f"{dir_out}/logs/bigwig/{{uniq_sample}}_bamCoverage.log"
    threads: 8
    shell:
        """
        bamCoverage -b {input.bam} \
            -o {output.bw} \
            --normalizeUsing BPM \
            --binSize 10 \
            -p {threads} \
            &> {log}
        """


# Rule 7.3: TSS Enrichment Calculation and Plotting
rule tss_enrichment:
    input:
        bw = f"{dir_out}/bigwig/{{uniq_sample}}.bw",
        tss = config.get("tss_bed", "hg38_tss.bed")
    output:
        matrix = f"{dir_out}/tss/{{uniq_sample}}_tss_matrix.gz",
        plot = f"{dir_out}/tss/{{uniq_sample}}_tss_enrichment.png",
        tab = f"{dir_out}/tss/{{uniq_sample}}_tss_profile.tab"
    log:
        f"{dir_out}/logs/tss/{{uniq_sample}}.log"
    threads: 8
    shell:
        """
        # 1. Calculate the signal matrix around TSS
        # -a: distance downstream of TSS
        # -b: distance upstream of TSS
        # --skipZeros: ignore regions with no signal
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 1000 -a 1000 \
            -R {input.tss} \
            -S {input.bw} \
            --skipZeros \
            -o {output.matrix} \
            -p {threads} &> {log}

        # 2. Generate the plot
        plotProfile \
            -m {output.matrix} \
            -o {output.plot} \
            --plotTitle "TSS Enrichment: {wildcards.uniq_sample}" \
            --regionsLabel "TSS" \
            --plotType lines \
            --outFileNameData {output.tab} \
            --perGroup \
            --colors green &>> {log}
        """


# # Rule 7.3.2: Prepare matrix for R
# rule tss_tab:
#     input:
#         matrix = f"{dir_out}/tss/{{uniq_sample}}_tss_matrix.gz"
#     output:
#  # Data for R
#     log:
#         f"{dir_out}/logs/tss/{{uniq_sample}}.tab.log"
#     threads: 8
#     shell:
#         """
#         # This converts the matrix into a simple profile table R can read
#         plotProfile -m {input.matrix} --outFileName {output.tab} &>> {log}
#         """

# Rule 7.3.3: Plot for all the samples
rule plot_all_tss:
    input:
        tabs = expand(f"{dir_out}/tss/{{uniq_sample}}_tss_profile.tab", uniq_sample=UNIQ_SAMPLES)
    output:
        plot = f"{dir_out}/plots/all_samples_tss_rich.png"
    log:
        f"{dir_out}/logs/plots/tss_rich_plot.log"
    params:
        tss_dir = f"{dir_out}/tss"
    shell:
        "Rscript plot_all_tss.R {output.plot} {params.tss_dir} &> {log}"
