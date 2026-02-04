#!/bin/bash
#SBATCH --job-name=spaf_spurio
#SBATCH --output=spaf_spurio.out
#SBATCH --error=spaf_spurio.err
#SBATCH --time=48:00:00
#SBATCH --mem=12GB

#rm output/summaries/query_summary.txt
python -u run_spurio.py part3_gpc_testing/sequences/spaf_sequences.fasta
#mv output/summaries/query_summary.txt output/summaries/spaf_summary.txt
cp spaf_sequences_summary.txt part3_gpc_testing/results/spaf/spaf_sequences_summary.txt
