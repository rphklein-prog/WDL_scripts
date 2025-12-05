version 1.0

# This workflow performs EdgeR differential analysis after salmon quantification, runs 
# gene ontology with Goseq and generates plots including mds, heatmap, and volcano plots.


workflow run_txi_edger {
  input {
    String transcriptome_name
    Array[String] data_groups
    Array[File] sf_samples
    String output_name
    File rscript
  }
  
  call run_edger {
    input:
      transcriptome_name = transcriptome_name,
      data_groups = data_groups,
      sf_samples = sf_samples,
      output_name = output_name,
      rscript = rscript
  }

  output {
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
  }
}
