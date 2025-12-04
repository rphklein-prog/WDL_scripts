#!/usr/bin/env Rscript

# this R script takes .sf files after salmon quantification of RNA-seq data, runs edgeR,
# and gene ontology analysis, and then creates various plots including mds, heatmap, and
# volcano plots. It takes the following arguments: first is transcriptome name 
# (ex EnsDb.Hsapiens.v75), the second is a file of the data groups (example: "Control", 
# "Control", "Treatment", "Treatment"), the third through second to last are salmon quant 
# file names, last is the name for the output file

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
stop("Usage: Rscript script.R <transcriptome_name> <data_groups_file> <salmon_files...> <output_file>")
}

transcriptome_name <- args[1]
DataGroups <- scan(args[2], what = "", sep = ",", quote = "\"", strip.white = TRUE)
files <- args[3:(length(args)-1)]
output_file <- args[length(args)]

if (length(DataGroups) != length(files)) {
    stop("Number of group labels does not match number of quant files.")
}

#import all necessary libraries

suppressPackageStartupMessages({
library(ensembldb)
library(GenomicFeatures)
library(transcriptome_name, character.only = TRUE)
library(tximport)
library(edgeR)
library(goseq)
library(GO.db)
library(org.Hs.eg.db)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ghibli)
})


# get names of samples and create table with a column for each sample's quant data for each
# transcript the transcript file
names(files) <- paste0("sample", seq_along(files))
edb <- get(transcriptome_name, asNamespace(transcriptome_name))
Tx <- transcripts(edb, return.type="DataFrame")
tx2gene <- subset(Tx, select=c("tx_name", "gene_id"))
txi.salmon <- tximport(files, type = "salmon", ignoreTxVersion=TRUE, tx2gene = tx2gene)

# create a DGEList object
group <- factor(DataGroups)
cts <- txi.salmon$counts
x <- DGEList(cts,group=group)
x.full <- x
design <- model.matrix(~ group)

# remove data with low counts in 2 or more samples
print(apply(x$counts, 2, sum))
keep <- rowSums(cpm(x)>5) >=2
x <- x[keep,]
# write.table(x$counts, file = paste0(output_file, ".counts.txt"), sep="\t", quote = FALSE, row.names = TRUE)
x$samples$lib.size <- colSums(x$counts)


#calculate normalization factors
x <-calcNormFactors(x)
x <- estimateDisp(x, design)
counts_table <- cpm(x, log=TRUE)
write.table(counts_table, file = paste0(output_file, ".logcounts.txt"), sep="\t", quote = FALSE, row.names = TRUE)


#run glm with quasi-likelihood
x1 <- glmQLFit(x, design)
x2 <- glmQLFTest(x1, coef=2)

#return top tags and write datatable
edgeR_table <- topTags(x2, n=Inf)$table
write.table(edgeR_table, file = paste0(output_file, ".txt"), sep = "\t", quote = FALSE, row.names = TRUE)
pdf(file= paste0(output_file, "_mds.pdf"))
plotMDS(x, method="bcv", col=as.numeric(x$samples$group))
dev.off()

genome_df <- data.frame(
  Ensembl= c("EnsDb.Hsapiens.v86", "EnsDb.Hsapiens.v75", "EnsDb.Mmusculus.v79"),
  Genome= c("hg38", "hg19", "mm10"))

## start GO-seq gene ontology analysis with data from EdgeR after qL calculation chose the
## significantly differentially expressed genes (up and down, just up, just down)and put 
##cthem in a new table

# convert transcript lengths → gene lengths (take median per gene)

for(subset in c("all", "upreg", "downreg")) {
  
  if(subset == "all") genes_vec <- edgeR_table$logFC != 0
  if(subset == "upreg") genes_vec <- edgeR_table$logFC > 0
  if(subset == "downreg") genes_vec <- edgeR_table$logFC < 0

  genes <- as.integer(p.adjust(edgeR_table$PValue[genes_vec], method="BH") < 0.05)
  names(genes) <- rownames(edgeR_table)[genes_vec]

  if(sum(genes) > 0) {
    pwf <- nullp(genes, "hg19", "ensGene")
    GO.wall <- goseq(pwf, "hg19", "ensGene") 
    # Adjust GO p-values
    GO.wall$padj <- p.adjust(GO.wall$over_represented_pvalue, method="BH")
    enriched.GO <- GO.wall$category[GO.wall$padj < 0.05]

    # Save output
    write.table(GO.wall, file = paste0(output_file, "_", subset, "_GO.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  }
}


# convert ensembl gene ids to gene symbol, if no gene symbol use ensembl gene id, also
# make sure each is unique 
edgeR_counts <- as.data.frame(counts_table)
symbols <- mapIds(org.Hs.eg.db, keys=rownames(edgeR_counts), keytype = "ENSEMBL", column = "SYMBOL")
symbols[is.na(symbols)] <- rownames(edgeR_counts)[is.na(symbols)]
symbols <- make.unique(symbols)
rownames(edgeR_counts) <- symbols

edgeR_fold <- as.data.frame(edgeR_table)
symbols <- mapIds(org.Hs.eg.db, keys=rownames(edgeR_fold), keytype = "ENSEMBL", column = "SYMBOL")
symbols[is.na(symbols)] <- rownames(edgeR_fold)[is.na(symbols)]
symbols <- make.unique(symbols)
rownames(edgeR_fold) <- symbols
# return only significant genes with FDR < 0.05
sig_genes <- edgeR_fold[edgeR_fold$FDR < 0.05, ]
write.table(sig_genes, file = paste0(output_file, "_sig_DEG.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# create a heatmap for the 100 most DEG
#add column names to new dataframe
colnames(edgeR_counts) <- DataGroups

# sort from highest logFC to lowest and select the top 100 rows
edgeR_sorted <-sig_genes[order(abs(sig_genes$logFC), decreasing=TRUE), ]
top100 <-edgeR_sorted[1:100, ]
subset_data <- edgeR_counts[rownames(edgeR_counts) %in% rownames(top100), ]

# plot heatmap
rnames <- rownames(subset_data)
mat_data <- data.matrix(subset_data)
rownames(mat_data) <- rnames
n_colors_per_segment <- 100
green_colors <- colorRampPalette(c("green", "black"))(n_colors_per_segment)
red_colors <- colorRampPalette(c("black", "red"))(n_colors_per_segment)
my_palette <- c(green_colors, red_colors)
col_breaks <- seq(min(mat_data), max(mat_data), length.out = length(my_palette) + 1)
pdf(file= paste0(output_file, "_heatmap.pdf"))
heatmap.2(mat_data, main = "Top Differentially Expressed Genes", notecol="NA", density.info="none", trace="none", col=my_palette, breaks=col_breaks, dendrogram="none", Colv=NA)
dev.off()

# make volcano plot
ghibli_subset <- ghibli_palette("SpiritedMedium")[c(1, 5)]

edgeR_fold$topDE <- "Neutral"
edgeR_fold$topDE[edgeR_fold$logFC > 0.25 & edgeR_fold$FDR < 0.05] <- "Up"
edgeR_fold$topDE[edgeR_fold$logFC < -0.25 & edgeR_fold$FDR < 0.05] <- "Down"

top_up <- edgeR_fold[edgeR_fold$topDE == "Up", ]
top_up <- top_up[order(top_up$FDR), ][1:min(5, nrow(top_up)), ]
top_down <- edgeR_fold[edgeR_fold$topDE == "Down", ]
top_down <- top_down[order(top_down$FDR), ][1:min(5, nrow(top_down)), ]
top_genes <- rbind(top_up, top_down)

pdf(file = paste0(output_file, "_volcano.pdf"))
ggplot(edgeR_fold, aes(x=logFC, y=-log10(FDR), color=topDE)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(
    data=top_genes,
    aes(label=row.names(top_genes)),
    size=4,
    max.overlaps = Inf
  ) +
  xlim(-3, 3) +
  ylim(0, max(-log10(edgeR_fold$FDR), na.rm=TRUE) + 1) +
  theme_bw() +
  theme(
    text = element_text(size=18),
    axis.line = element_line(color="darkgrey", size=1),
    axis.ticks = element_line(colour="black", size=1)
  ) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Neutral" = "grey")) +
  labs(
    title="Volcano plot of Differential Expression",
    x="Log2 Fold Change",
    y="-log10(FDR)",
    color="DE Status"
  )
dev.off()
