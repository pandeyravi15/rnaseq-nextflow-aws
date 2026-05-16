#!/usr/bin/env Rscript

# Get arguments
args <- commandArgs(trailingOnly=TRUE)
counts_file  <- args[1]
condition_ad <- args[2]
pval_cutoff  <- as.numeric(args[3])
lfc_cutoff   <- as.numeric(args[4])

suppressMessages({
    library(DESeq2)
    library(org.Mm.eg.db)
    library(ggplot2)
})

cat("=== DESeq2 Analysis Starting ===\n")

# Load counts
counts <- read.table(
    counts_file,
    header=TRUE, row.names=1, sep=""
)
cat("Genes:", nrow(counts), "| Samples:", ncol(counts), "\n")

# Create metadata
sample_names <- colnames(counts)
condition <- ifelse(grepl(condition_ad, sample_names), "AD", "WT")
coldata <- data.frame(
    row.names = sample_names,
    condition = factor(condition, levels = c("WT", "AD"))
)
cat("AD:", sum(condition=="AD"), "| WT:", sum(condition=="WT"), "\n")

# DESeq2
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = coldata,
    design    = ~ condition
)
keep <- rowSums(counts(dds) >= 10) >= 5
dds  <- dds[keep,]
cat("Genes after filtering:", nrow(dds), "\n")

cat("Running DESeq2...\n")
dds <- DESeq(dds)
res <- results(dds,
               contrast = c("condition","AD","WT"),
               alpha    = pval_cutoff)

# Add gene symbols
gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys      = rownames(res),
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
)
res_df            <- as.data.frame(res)
res_df$ensembl_id  <- rownames(res_df)
res_df$gene_symbol <- gene_symbols
res_df            <- res_df[order(res_df$padj),]

# Significant genes
sig <- res_df[!is.na(res_df$padj) & res_df$padj < pval_cutoff,]

cat("\n=== RESULTS ===\n")
cat("Significant genes:", nrow(sig), "\n")
cat("Upregulated:", sum(sig$log2FoldChange > 0, na.rm=TRUE), "\n")
cat("Downregulated:", sum(sig$log2FoldChange < 0, na.rm=TRUE), "\n")

# Save results
dir.create("deseq2_results", showWarnings=FALSE)
write.csv(res_df, "deseq2_results/all_genes.csv",         row.names=FALSE)
write.csv(sig,    "deseq2_results/significant_genes.csv", row.names=FALSE)
cat("Results saved!\n")

# Volcano plot
vd <- res_df[!is.na(res_df$padj),]
vd$significance <- "Not significant"
vd$significance[vd$padj < pval_cutoff &
                vd$log2FoldChange >  lfc_cutoff] <- "Upregulated"
vd$significance[vd$padj < pval_cutoff &
                vd$log2FoldChange < -lfc_cutoff] <- "Downregulated"

ad_genes  <- c("Trem2","Apoe","Bin1","Camk2a","Spp1",
               "Lpl","Tyrobp","Cst7","P2ry12","Cx3cr1")
vd$label <- ifelse(vd$gene_symbol %in% ad_genes,
                   vd$gene_symbol, NA)

p <- ggplot(vd, aes(x=log2FoldChange, y=-log10(padj),
                    color=significance, label=label)) +
    geom_point(alpha=0.6, size=1.5) +
    scale_color_manual(values=c(
        "Upregulated"     = "red",
        "Downregulated"   = "blue",
        "Not significant" = "grey60"
    )) +
    geom_text(hjust=-0.1, vjust=0, size=3,
              color="black", na.rm=TRUE) +
    geom_vline(xintercept=c(-lfc_cutoff, lfc_cutoff),
               linetype="dashed") +
    geom_hline(yintercept=-log10(pval_cutoff),
               linetype="dashed") +
    theme_classic() +
    labs(title="3xTgAD vs WT — Cortex DESeq2",
         x="Log2 Fold Change", y="-log10(padj)",
         color="Direction") +
    theme(plot.title=element_text(hjust=0.5, face="bold"))

ggsave("deseq2_results/volcano_plot.png",
       plot=p, width=10, height=8, dpi=300)

cat("Volcano plot saved!\n")
cat("=== Analysis Complete ===\n")
