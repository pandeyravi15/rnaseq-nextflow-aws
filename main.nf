#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================
// PROCESS 1 — Download counts matrix from S3
// ================================================
process DOWNLOAD_COUNTS {
    tag "Downloading counts matrix"

    output:
    path "counts_matrix.txt", emit: counts

    script:
    """
    echo "Downloading counts matrix from S3..."
    aws s3 cp s3://${params.bucket}/${params.counts_file} counts_matrix.txt
    echo "Done: \$(wc -l < counts_matrix.txt) lines"
    """
}

// ================================================
// PROCESS 2 — Run DESeq2 via R script
// ================================================
process RUN_DESEQ2 {
    tag "DESeq2 analysis"
    publishDir "${System.getenv("HOME")}/projects/rnaseq/results/nextflow_deseq2", mode: "copy"

    input:
    path counts_file

    output:
    path "deseq2_results/*.csv", emit: results
    path "deseq2_results/*.png", emit: plots

    script:
    """
    Rscript ${projectDir}/bin/deseq2_analysis.R \
        ${counts_file} \
        ${params.condition_ad} \
        ${params.pval_cutoff} \
        ${params.lfc_cutoff}
    """
}

// ================================================
// PROCESS 3 — Upload to S3
// ================================================
process UPLOAD_RESULTS {
    tag "Uploading to S3"

    input:
    path results
    path plots

    script:
    """
    echo "Uploading to S3..."
    aws s3 sync . s3://${params.bucket}/${params.outdir}/
    echo "Upload complete!"
    aws s3 ls s3://${params.bucket}/${params.outdir}/ --human-readable
    """
}

// ================================================
// WORKFLOW
// ================================================
workflow {
    log.info "Starting RNA-seq DESeq2 Pipeline — Ravi Pandey"

    DOWNLOAD_COUNTS()
    RUN_DESEQ2(DOWNLOAD_COUNTS.out.counts)
    UPLOAD_RESULTS(
        RUN_DESEQ2.out.results,
        RUN_DESEQ2.out.plots
    )

    log.info "Pipeline submitted!"
}
