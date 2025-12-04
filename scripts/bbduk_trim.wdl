version 1.0

#trims Illumina adapters from the ends of paired reads using bbmap program bbduk

workflow bbduk_trimming {
  input {
    File fastq1
    File fastq2
  }

  call bbduk {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2

  }

  output {
    File trimmed_read1 = bbduk.trimmed_read1
    File trimmed_read2 = bbduk.trimmed_read2
  }
}
# change adapter sequences in the bbduk task if trimming files from other sequencing 
# platforms

task bbduk {
  input {
    File fastq1
    File fastq2
    String adapter_F = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
    String adapter_R = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG"
  }
  
  #dynamically estimate needed memory
  
  Int disk_size = ceil((size(fastq1, "GB") + size(fastq2, "GB")) * 3)

  # Create an output directory, get basename from input file, can be .fastq, .fq, fq.gz,
  # or .fastq.gz
  
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
    File trimmed_read1 = "~{prefix}_trimmed_1.fastq.gz"
    File trimmed_read2 = "~{prefix}_trimmed_2.fastq.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/bbmap:38.96--h5c4e2a8_0"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
  }
}
