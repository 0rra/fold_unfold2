#!/bin/bash
#SBATCH --job-name=afdb_seqs
#SBATCH --time=24:00:00
#SBATCH --mem=20G

wget https://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta

mv sequences.fasta data/afdb_paes/afdb_sequences.fasta 