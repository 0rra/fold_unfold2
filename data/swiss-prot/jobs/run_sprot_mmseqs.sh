#!/bin/bash
#SBATCH --job-name=run_sprot_mmseqs
#SBATCH --output=logs/run_sprot_mmseqs.log
#SBATCH --error=logs/run_sprot_mmseqs.err
#SBATCH --time=36:00:00
#SBATCH --mem=50GB

source /homes/${USER}/.bashrc

sp_data_dir="data/swiss-prot/"

mmseqs easy-cluster ${sp_data_dir}/2025_03_uniprot_sprot.fasta ${sp_data_dir}/2025_03_uniprot_sprot_30 tmp_sprot --min-seq-id 0.3

rm -r tmp_sprot
