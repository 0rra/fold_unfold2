#!/bin/bash
#SBATCH --job-name=traf_1000_shadow
#SBATCH --output=traf_1000_shadow.out
#SBATCH --error=traf_1000_shadow.err
#SBATCH --time=45:00:00
#SBATCH --mem=12GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ako@ebi.ac.uk

source /homes/ako/.bashrc
conda activate pdb2
cd /nfs/research/agb/research/ako/git_repos/shadow_test
python -u /nfs/research/agb/research/ako/git_repos/shadow/shadow_test.py /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/sequences/traf_1000_sequences.fasta --resultsdir /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/results/traf_1000/shadow_res --mode diamond --csv --trembl --bacteria
