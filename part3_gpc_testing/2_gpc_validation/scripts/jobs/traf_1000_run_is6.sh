#!/bin/bash
#SBATCH --job-name=traf_1000_is6
#SBATCH --output=traf_1000_is6.log
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user='ako@ebi.ac.uk'

source /homes/${USER}/.bashrc
nextflow run ebi-pf-team/interproscan6 -r 6.0.0 	-profile singularity --datadir data 	--input /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/sequences/traf_1000_sequences.fasta --applications Pfam,AntiFam,mobidblite --formats tsv 	--outprefix traf_1000_is6
mv traf_1000_is6.tsv /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/results/traf_1000/traf_1000_is6.tsv
