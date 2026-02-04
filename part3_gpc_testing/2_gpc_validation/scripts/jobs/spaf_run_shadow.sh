#!/bin/bash
#SBATCH --job-name=spaf_shadow
#SBATCH --output=spaf_shadow.out
#SBATCH --error=spaf_shadow.err
#SBATCH --time=45:00:00
#SBATCH --mem=12GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ako@ebi.ac.uk

source /homes/${USER}/.bashrc
conda activate pdb2
cd /nfs/research/agb/research/ako/git_repos/shadow_test
python -u /nfs/research/agb/research/ako/git_repos/shadow/shadow_test.py /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/sequences/spaf_sequences.fasta --resultsdir /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/results/spaf/shadow_res --mode diamond --csv --trembl --bacteria
