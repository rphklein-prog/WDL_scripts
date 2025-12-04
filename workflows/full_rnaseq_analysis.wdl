version 1.0

# This workflow starts with a list of SRA accession numbers (can be a Terra data table with 
# SRR ids in a column), a list matching each SRA accession number to a sample type (ex- 
# Treatment or Control), and the desired salmon index for quantification. It pulls data 
# from sra then trims adapters with bbduk (defaults are Illumina universal adapters), runs
# fastqc to check that adapters have been removed and sequence quality is acceptable, 
# quantifies reads with salmon and then performs EdgeR differential analysis and generates
# plots. Assumes data is paired-end with F and R reads.


workflow full_rna_seq {
  meta {
    author: "Rachel Klein"
    description: "Workflow for analyzing RNA-seq data with Salmon and EdgeR"
    outputs: {
      quant_sfs: "array of salmon quantification files for all samples"
      log_counts: "table containing all samples with log normalized counts per million for each gene"
      edger_table: "output from edgeR after model fitting and differential analysis"
      mds_plot: "multi-dimensional scaling plot of analyzed samples"
      go_all: "output of GOseq analysis for all differentially expressed genes"
      go_up: "output of GOseq analysis for upregulated genes"
      go_down: "output of GOseq analysis for downregulated genes"
      sig_genes: "file of significantly differentially expressed genes (subset of edger_table)"
      heatmap: "heatmap of 100 most differentially expressed genes (by logFC)"
      volcano: "volcano plot of edgeR results"
    }
  }
  
  parameter_meta {
    sra_accesions: "list of SRA sample ID's to be pulled and analyzed"
    groups: "list IDing the group for each sample (ex Control, Control, Treatment, ...)"
    addldisk: "increase disk size by this much for running fastqc"
    cpu: "number of CPUs to use for fastqc"
    memory: "memory allocation in GB for fastqc"
    salmon_index_tar: "tar compressed file of salmon index, can be downloaded from refgenie and added to terra strorage bucket for referencing in workflow"
    transcriptome_name: "transcriptome to use for tximports to generate counts table with ids"
    output_name: "name for the output files and plots"
    rscript: "Full_R_analysis.r rscript to run edgeR and plotting functions"
  }
  
  
    
  input {
    Array[String] sra_accessions
    Array[String] groups
    Int addldisk = 10
	Int cpu = 4
	Int memory = 8
	File salmon_index_tar
	String transcriptome_name
    String output_name
    File rscript= "https://raw.githubusercontent.com/rphklein-prog/WDL_scripts/main/Full_R_analysis.r" 
  }
  
  scatter (sra in sra_accessions){
    call fastq_dl_sra {
      input:
        sra_accession=sra
    }

    # use bbduk to trim adapters from reads (default are Illumina universal adapters)
    call bbduk {
      input:
        fastq1 = fastq_dl_sra.reads[0],
        fastq2 = fastq_dl_sra.reads[1] 
    }
  
    # run fastqc on both forward and reverse reads by scattering
    scatter (fq in bbduk.trimmed_reads){
      call fastqc {
   	    input:
          fastq = fq,
          addldisk = addldisk,
          cpu = cpu,
          memory = memory    
      }
    }
  
    call salmon_quant {
      input:
        fastq1 = bbduk.trimmed_reads[0],
        fastq2 = bbduk.trimmed_reads[1],
        sra = sra,
        index_tar = salmon_index_tar
    }
}

# Flatten nested outputs
  Array[Array[Array[File]]] nested_fastqc = fastqc.html_report
  Array[Array[File]] mid_fastqc = flatten(nested_fastqc)
  Array[File] fastqc_reports = flatten(mid_fastqc)
  
  Array[File] quant_sfs_all = salmon_quant.quant_sf
  
  call run_edger {
    input:
      transcriptome_name = transcriptome_name,
      data_groups = groups,
      sf_samples = quant_sfs_all,
      output_name = output_name,
      rscript = rscript
  }  
  
  output {
    Array[File] all_fastqc_reports = fastqc_reports
    Array[File] quant_sfs = quant_sfs_all
    File log_counts = run_edger.log_counts
    File edger_table = run_edger.edger_table 
    File mds_plot = run_edger.mds_plot
    File go_all = run_edger.go_all
    File go_up = run_edger.go_up
    File go_down = run_edger.go_down
    File sig_genes = run_edger.sig_genes
    File heatmap = run_edger.heatmap
    File volcano = run_edger.volcano
  }
}

task fastq_dl_sra {
  input {
    String sra_accession
  }
  
  command <<<
    # write version to VERSION file for records
    fastq-dl --version | tee VERSION
    fastq-dl --accession ~{sra_accession}

    # tag single-end reads with _1
    if [ -f "~{sra_accession}.fastq.gz" ] && [ ! -f "~{sra_accession}_1.fastq.gz" ]; then
      mv "~{sra_accession}.fastq.gz" "~{sra_accession}_1.fastq.gz"
    fi
  >>>
  
  output {
    Array[File] reads = glob("~{sra_accession}_*.fastq.gz")
  }
  
  runtime {
    docker: "quay.io/biocontainers/fastq-dl:3.0.1--pyhdfd78af_0"
    memory:"8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible:  1
  }
}

task fastqc {
	input {
		File fastq
		Int addldisk = 10
		Int cpu = 4
		Int memory = 8
	}
	# dynamically calculate disk needs
	Int finalDiskSize = addldisk + ceil(size(fastq, "GB"))

	command <<<
		mkdir outputs
		fastqc -o outputs ~{fastq}
	>>>

	runtime {
		cpu: cpu
		docker: "biocontainers/fastqc:v0.11.9_cv8"
		disks: "local-disk " + finalDiskSize + " SSD"
		memory: memory + " GB"
		preemptible: 1
	}

	output {
		Array[File] html_report = glob("outputs/*_fastqc.html")
	}
	
}

task bbduk {
  input {
    File fastq1
    File fastq2
    String adapter_F = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    String adapter_R = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG"
  }
  
  # dynamically calculate disk needs
  Int disk_size = ceil((size(fastq1, "GB") + (if defined(fastq2) then size(fastq2, "GB") else 0)) * 3)
  
  # strip the name from the fastq1 file to use for outputs
  String prefix = sub(basename(fastq1), "\\.f(ast)?q(\\.gz)?$", "")
  
  
  command <<<
    set -euxo pipefail
    # Run bbduk
    bbduk.sh \
      in1=~{fastq1} \
      in2=~{fastq2} \
      out1=~{prefix}_trimmed_1.fastq.gz \
      out2=~{prefix}_trimmed_2.fastq.gz \
      literal=~{adapter_F},~{adapter_R} \
      ktrim=r \
      k=23 \
      mink=11 \
      hdist=1 \
      tpe \
      tbo
  >>>

  output {
    # Match files for both forward and reverse reads if present
    Array[File] trimmed_reads = glob("~{prefix}_trimmed_*.fastq.gz")
  }

  runtime {
    docker: "quay.io/biocontainers/bbmap:38.96--h5c4e2a8_0"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
  }
}


# Run salmon quantification
task salmon_quant {
  input {
    File fastq1
    File fastq2
    String sra
    File index_tar
  }

  command <<<
    set -euxo pipefail
    
    # Extract index locally
    tar -xzf ~{index_tar}
    INDEX_DIR=$(tar -tzf ~{index_tar} | head -1 | cut -f1 -d"/" || echo ".")
    
    salmon quant \
      -i "${INDEX_DIR}" \
      -l A \
      -1 ~{fastq1} \
      -2 ~{fastq2} \
      -p 4 \
      -o "~{sra}_quant" \
      --seqBias \
      --useVBOpt \
      --validateMappings 
  >>>
		
  output {
    File quant_sf = "~{sra}_quant/quant.sf"
    Array[File] quant_results = glob("~{sra}_quant/*")
  }
	
  runtime {
    docker: "quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"
    cpu: 4
    memory: "32G"
    disks: "local-disk 200 SSD"  # Increased for large index
    preemptible: 1
  }
}

task run_edger {
  input {
    String transcriptome_name
    Array[String] data_groups
    Array[File] sf_samples
    String output_name
    File rscript
  }
  
  command <<<
    set -euo pipefail
    DATA_GROUPS_FILE="data_groups.txt"
    printf "%s\n" ~{sep="," data_groups} > $DATA_GROUPS_FILE
    SALMON_FILES="~{sep=' ' sf_samples}"
    Rscript ~{rscript} ~{transcriptome_name} $DATA_GROUPS_FILE $SALMON_FILES ~{output_name}
    >>>

  output {
    File log_counts = "~{output_name}.logcounts.txt"
    File edger_table = "~{output_name}.txt"
    File mds_plot = "~{output_name}_mds.pdf"
    File go_all = "~{output_name}_all_GO.txt"
    File go_up = "~{output_name}_upreg_GO.txt"
    File go_down = "~{output_name}_downreg_GO.txt"
    File sig_genes = "~{output_name}_sig_DEG.txt"
    File heatmap = "~{output_name}_heatmap.pdf"
    File volcano = "~{output_name}_volcano.pdf"

  } 

runtime {
    docker: "rphklein/txi_edger_plots:1.0"
    cpu: 2
    memory: "4G"
    preemptible: 1
  }
}
