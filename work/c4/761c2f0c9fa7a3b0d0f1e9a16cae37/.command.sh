#!/bin/bash -ue
Rscript /home/ec2-user/projects/nextflow_pipeline/bin/deseq2_analysis.R         counts_matrix.txt         3xTgAD         0.05         0.5
