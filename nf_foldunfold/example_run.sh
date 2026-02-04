#!/bin/bash
#SBATCH --job-name=af_pred
#SBATCH --output=af_pred_%a.log
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:a100:1
#SBATCH --array=1-4
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=#EMAIL HERE

outputdir="structure_results/antifam_structures"
inputdir="input_sequences/antifam_seed_seqs"
i=${SLURM_ARRAY_TASK_ID}

module load nextflow

nextflow run nf_foldunfold/main.nf \
--predictmethods alphafold3,colabfold,esmfold \
--input ${inputdir}/seed_seqs_part_00${i}.fasta \
--outdir ${outputdir}
