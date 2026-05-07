# CUT&Tag-seq preprocessing with snakemake 

* Initial QC with fastqc and multiqc
* Trimming with fastp
* Post-trimming QC with fastqc and multiqc
* Alignment to hg38 with bowtie2
* Remove mitochondrial reads (and plot)
* Remove ENCODE blacklist regions
* Remove duplicate reads with samtools markdup
* Peak calling with macs3
* FRiP scores and plot
* Calculate TSS enrichment and plot
* Generate bigWig tracks for visualization

