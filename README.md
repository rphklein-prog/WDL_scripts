# WDL_scripts
WDL scripts for bioinformatic analysis on Terra. Included scripts:

fetch_trim_fastqc_salmon.wdl: combines functions into a single workflow: fetches data from SRA repository based on SRR accesions, 
trims adapters from reads, runs fastqc, and then quantifies counts with salmon.

txi_edger_go_plot.wdl: a WDL wrapper to run provided Rscript (Full_analysis.R). Takes salmon output (quant.sf files), uses tximports
to generate a counts table, runs edgeR, performs GO analysis, and generates plots (heatmap, mds, and volcano).

bbduk_trim: runs bbduk (from bbmap) to trim given adapter sequences from the ends of reads. See https://github.com/BioInfoTools/BBMap

salmon_wdl: runs salmon read quantification on paired-end RNA-seq data. See https://salmon.readthedocs.io/en/latest/salmon.html




