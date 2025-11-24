version 1.0

# This workflow starts with a Terra data table with SRR ids in a column. It pulls data 
# from sra then trims adapters with bbduk (defaults are Illumina universal adapters), runs
# fastqc to check that adapters have been removed and sequence quality is acceptable, and 
# finally, quantifies reads with salmon. 

workflow fetch_sra_to_fastqc {
  input {
    String sra_accession
    Int addldisk = 10
	Int cpu = 4
	Int memory = 8
	String sample_id
	File salmon_index_tar
  }
  
  call fastq_dl_sra {
    input:
      sra_accession=sra_accession
  }

# use bbduk to trim adapters from reads (default are Illumina universal adapters)
  call bbduk {
    input:
      fastq1 = fastq_dl_sra.reads[0]
      fastq2 = if size(fastq_dl_sra.reads) > 1 then fastq_dl_sra.reads[1] else null  
  }
  
  # run fastqc on both forward and reverse reads by scattering
  scatter (fq in bbduk.trimmed_reads){
    call fastqc {
   	  input:
        fastq = fq
        addldisk = addldisk
        cpu = cpu
        memory = memory    
    }
  }
    
  call salmon_quant {
    input:
      fastq1 = bbduk.trimmed_reads[0]
      fastq2 = if size(bbduk.trimmed_reads) > 1 then bbduk.trimmed_reads[1] else null
      salmon_index_tar = salmon_index_tar,
      sample_id = sample_id
  }

  
  output {
    Array[File] reads = fastq_dl_sra.reads
    Array[File] report = flatten(fastqc.*.html_report)
    Array[File] trimmed_reads = bbduk.trimmed_reads
    File quant_sf = salmon_quant.quant_sf
    Array[File] quant_results = salmon_quant.quant_results
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
    File? fastq2
    String adapter_F = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    String adapter_R = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG"
  }
  
  # dynamically calculate disk needs
  Int disk_size = ceil((size(fastq1, "GB") + (fastq2 ? size(fastq2, "GB") : 0)) * 3)
  
  # strip the name from the fastq1 file to use for outputs
  String prefix = sub(basename(fastq1), "\\.f(ast)?q(\\.gz)?$", "")
  
  # set up to handle situation with only 1 fastq (single end sequencing)
  String out1 = "~{prefix}_trimmed_1.fastq.gz"
  String? out2 = fastq2 ? "~{prefix}_trimmed_2.fastq.gz"
  
  command <<<
    set -euxo pipefail
    # Run bbduk
    if [ -f "~{fastq2}" ]; then
      bbduk.sh \
      in1=~{fastq1} \
      in2=~{fastq2} \
      out1=~{out1} \
      out2=~{out2} \
      literal=~{adapter_F},~{adapter_R} \
      ktrim=r \
      k=23 \
      mink=11 \
      hdist=1 \
      tpe \
      tbo
    else
      bbduk.sh \
      in=~{fastq1} \
      out=~{out1} \
      literal=~{adapter_F},~{adapter_R} \
      ktrim=r \
      k=23 \
      mink=11 \
      hdist=1
    fi
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
    File? fastq2
    String sample_id
    File salmon_index_tar
  }

# run salmon	
  command <<<
    set -euo pipefail
    # extract salmon index and put into a folder
    tar -xzvf ~{salmon_index_tar}
    INDEX_DIR=$(tar -tzf ~{salmon_index_tar} | head -1 | cut -f1 -d"/" || echo ".")
    # use reverse reads if present, otherwise quantify as single end
    if [ -f "~{fastq2}" ]; then
      salmon quant \
        -i "${INDEX_DIR}" \
        -l A \
        -1 ~{fastq1} \
        -2 ~{fastq2} \
        -p 4 \
        -o "~{sample_id}_quant" \
        --seqBias \
        --useVBOpt \
        --validateMappings 
    else
      salmon quant \
        -i "${INDEX_DIR}" \
        -l A \
        -r ~{fastq1} \
        -p 4 \
        -o "~{sample_id}_quant" \
        --seqBias \
        --useVBOpt \
        --validateMappings 
    fi
  >>>
		
  output {
    File quant_sf = "~{sample_id}_quant/quant.sf"
    Array[File] quant_results = glob("~{sample_id}_quant/*")
  }

# default runtime set-up, modify as needed		
  runtime {
    docker: "quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5"
    cpu: 4
    memory: "32G"
    disks: "local-disk 250 SSD"
    preemptible: 1
  }
}
