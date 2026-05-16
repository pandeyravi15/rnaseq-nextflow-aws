#!/bin/bash -ue
echo "Downloading counts matrix from S3..."
aws s3 cp s3://ravi-pandey-bioinformatics/rnaseq/results/counts/5XFAD_counts.txt counts_matrix.txt
echo "Done: $(wc -l < counts_matrix.txt) lines"
