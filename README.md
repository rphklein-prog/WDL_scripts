# WDL_scripts
WDL workflows and scripts for RNA-seq analysis on Terra. Included scripts:

full_rnaseq_analysis.wdl: combines functions into a single workflow: starts with a list of SRA accession numbers (can be a Terra data table with 
SRR ids in a column), a list matching each SRA accession number to a sample type (ex- Treatment or Control), and the desired salmon index for 
quantification. It pulls data from sra then trims adapters with bbduk (defaults are Illumina universal adapters), runs fastqc to check that 
adapters have been removed and sequence quality is acceptable, quantifies reads with salmon and then performs EdgeR differential analysis and 
generates plots. Assumes data are paired-end with F and R reads.

fetch_trim_fastqc_salmon.wdl: fetches data from SRA repository based on SRR accessions, trims adapters from reads, runs fastqc, and then 
quantifies counts with salmon. The expected input is a terra table with an id column of SRA accession numbers and a path to the data bucket with 
the salmon index.

txi_edger_go_plot.wdl: a WDL wrapper to run provided Rscript (Full_R_analysis.r). Takes salmon output (quant.sf files), uses tximports
to generate a counts table, runs edgeR, performs GO analysis, and generates plots (heatmap, mds, and volcano). Requires as input the
quant.sf files for the samples (ideally in one column of a terra data table), the name of the reference transcriptome (ie- 
EnsDb.Hsapiens.v75), a text file or list of the sample types in order ("Control", "Control", "Treatment", "Treatment"), the output name
for the files, and the Rscript location.

bbduk_trim.wdl: runs bbduk (from bbmap) to trim given adapter sequences from the ends of reads. See https://github.com/BioInfoTools/BBMap

fastqcWF.wdl: runs fastqc on files and provides stats on read quality. See https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

salmon_wdl.wdl: runs salmon read quantification on paired-end RNA-seq data. See https://salmon.readthedocs.io/en/latest/salmon.html

From R script:

Tximport: generate a counts table from quantification files (https://bioconductor.org/packages/release/bioc/html/tximport.html)

EdgeR: differential expression analysis (https://bioconductor.org/packages/release/bioc/html/edgeR.html)

GOseq: perform gene ontology analysis (https://bioconductor.org/packages/release/bioc/html/goseq.html)

Also included for reference- Dockerfile for Docker image rphklein/txi_edger_plots:1.0 used for txi_edger_go_plot.wdl (image can 
be pulled from Dockerhub)
