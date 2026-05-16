#!/bin/bash -ue
echo "Uploading to S3..."
aws s3 sync . s3://ravi-pandey-bioinformatics/rnaseq/results/nextflow/
echo "Upload complete!"
aws s3 ls s3://ravi-pandey-bioinformatics/rnaseq/results/nextflow/ --human-readable
