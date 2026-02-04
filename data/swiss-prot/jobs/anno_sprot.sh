#!/bin/bash
#SBATCH --job-name=anno_sprot2
#SBATCH --output=logs/anno_sprot2_%A.log
#SBATCH --error=logs/anno_sprot2_%A.err
#SBATCH --time=36:00:00
#SBATCH --mem=10GB

source /homes/${USER}/.bashrc
conda activate foldunfold2

fasta_dir="data/swiss-prot"

python -u data/swiss-prot/scripts/anno_sprot.py ${fasta_dir}/2025_03_uniprot_sprot.fasta ${fasta_dir}/2025_03_uniprot_sprot

# (optional) update sequence headers to only include ids (header descriptions may not be compatible with other tools/ result in parsing issues)
python data/swiss-prot/scripts/fix_names.py ${fasta_dir}/2025_03_uniprot_sprot.fasta > ${fasta_dir}/2025_03_uniprot_sprot_head_fixed.fasta
