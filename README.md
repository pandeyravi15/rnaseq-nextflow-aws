# RNA-seq DESeq2 Nextflow Pipeline on AWS

End-to-end automated RNA-seq differential expression pipeline built with Nextflow DSL2, 
running on AWS EC2/Batch with S3 for data storage.

## Overview

This pipeline performs differential expression analysis of RNA-seq count data using DESeq2, 
with full automation from data download through results upload to S3.

**Biological application:** Alzheimer's disease mouse model transcriptomics  
**Dataset:** 3xTgAD vs WT cortex (GEO: GSE161904)  
**Author:** Ravi Shanker Pandey, Ph.D. | The Jackson Laboratory

## Pipeline Architecture
S3 (input counts matrix)
↓
DOWNLOAD_COUNTS  — pulls data from AWS S3
↓
RUN_DESEQ2       — differential expression analysis
↓
UPLOAD_RESULTS   — pushes results back to S3

## Key Features

- **Nextflow DSL2** — modular, scalable pipeline architecture
- **AWS integration** — EC2 compute, S3 storage, Batch-ready
- **Automated DESeq2** — complete differential expression workflow
- **Gene annotation** — ENSEMBL to gene symbol mapping via org.Mm.eg.db
- **Publication-quality plots** — volcano plots with AD-relevant gene labels
- **Reproducible** — Docker/Singularity compatible, timestamped logging

## Pipeline Components

| File | Description |
|---|---|
| `main.nf` | Main Nextflow DSL2 workflow |
| `nextflow.config` | Pipeline configuration (local + AWS Batch profiles) |
| `bin/deseq2_analysis.R` | DESeq2 differential expression R script |

## Requirements

- Nextflow >= 23.0
- Java >= 17
- R >= 4.3 with DESeq2, ggplot2, org.Mm.eg.db
- AWS CLI configured with S3 access
- Docker (for containerized execution)

## Usage

### Local execution
```bash
nextflow run main.nf -profile local
```

### AWS Batch execution
```bash
nextflow run main.nf -profile aws
```

### Custom parameters
```bash
nextflow run main.nf -profile local \
  --bucket your-s3-bucket \
  --counts_file path/to/counts.txt \
  --pval_cutoff 0.05 \
  --lfc_cutoff 0.5
```

## Results

The pipeline generates:
- `all_genes.csv` — DESeq2 results for all tested genes
- `significant_genes.csv` — Filtered significant genes (padj < 0.05)
- `volcano_plot.png` — Publication-quality volcano plot

## Biological Results

Applied to 3xTgAD vs WT cortex data (30 samples, 43,629 genes):
- **2,339 significant genes** identified (padj < 0.05)
- **Upregulated in AD:** Spp1, Lpl, Tyrobp, Hexb (neuroinflammation)
- **Downregulated in AD:** Bin1, Camk2a (synaptic dysfunction)

## Author

**Ravi Shanker Pandey, Ph.D.**  
Associate Computational Scientist, The Jackson Laboratory  
[ORCID](https://orcid.org/0000-0001-9567-2851) | [GitHub](https://github.com/pandeyravi15) | [LinkedIn](https://linkedin.com/in/ravi-pandey-50943980)
