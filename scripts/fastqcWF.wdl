version 1.0
# runs fastqc

workflow FastqcWF {
	input {
		File fastq
	}

	call fastqc {
		input:
			fastq = fastq,
	}

	output {
		File report = fastqc.html_report
	}

}

task fastqc {
	input {
		File fastq
		Int addldisk = 10
		Int cpu = 4
		Int memory = 8
		Int preempt = 1 
	}
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
		preemptible: preempt
	}

	output {
		File html_report = glob("outputs/*_fastqc.html")[0]
	}
	
}
