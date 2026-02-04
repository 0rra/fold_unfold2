#!/bin/bash
#SBATCH --job-name=run_sprot_is6
#SBATCH --output=logs/run_sprot_is6.log
#SBATCH --error=logs/run_sprot_is6.err
#SBATCH --time=20:00:00
#SBATCH --mem=100GB

source /homes/${USER}/.bashrc

fasta_dir="data/swiss-prot"
module load nextflow/25.04.6
nextflow run ebi-pf-team/interproscan6 -r 6.0.0	-profile singularity --datadir data 	\\
--input ${fasta_dir}/2025_03_uniprot_sprot.fasta --applications Pfam,AntiFam,MobiDB-lite --formats tsv 	--outprefix 2025_03_uniprot_sprot_is6
