#!/bin/bash
#SBATCH --job-name=anno_sprot_tools
#SBATCH --output=logs/anno_sprot_tools.log
#SBATCH --error=logs/anno_sprot_tools.err
#SBATCH --time=36:00:00
#SBATCH --mem=50GB

source /homes/${USER}/.bashrc
conda activate foldunfold2

sp_data_dir="data/swiss-prot"

python parse_is6.py ${sp_data_dir}/2025_03_uniprot_sprot_is6.tsv ${sp_data_dir}/2025_03_uniprot_sprot_anno.tsv ${sp_data_dir}/2025_03_uniprot_sprot_anno.tsv
echo 'parsed is6'
python anno_clusters.py ${sp_data_dir}/2025_03_uniprot_sprot_30_cluster.tsv ${sp_data_dir}/2025_03_uniprot_sprot_anno.tsv ${sp_data_dir}/2025_03_uniprot_sprot_anno.tsv 
echo 'parsed mmseqs clusters'