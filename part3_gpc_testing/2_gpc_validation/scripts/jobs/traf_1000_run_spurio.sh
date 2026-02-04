#!/bin/bash
#SBATCH --job-name=traf_1000_spurio
#SBATCH --output=traf_1000_spurio.out
#SBATCH --error=traf_1000_spurio.err
#SBATCH --time=48:00:00
#SBATCH --mem=12GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ako@ebi.ac.uk

source /homes/ako/.bashrc
source /homes/ako/.bashrc
conda activate /nfs/research/agb/research/ako/software/miniforge3/envs/spurio_env2
cd /hps/software/users/agb/research/ako/spurio
#mv output/summaries/query_summary.txt
python -u run_spurio.py /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/sequences/traf_1000_sequences.fasta
#mv output/summaries/query_summary.txt output/summaries/traf_1000_summary.txt
cp traf_1000_sequences_summary.txt /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/results/traf_1000/traf_1000_sequences_summary.txt
