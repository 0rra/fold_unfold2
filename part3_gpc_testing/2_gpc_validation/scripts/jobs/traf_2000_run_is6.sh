#!/bin/bash
#SBATCH --job-name=traf_2000_is6
#SBATCH --output=traf_2000_is6.log
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user='ako@ebi.ac.uk'

source /homes/ako/.bashrc
nextflow run ebi-pf-team/interproscan6 -r 6.0.0 	-profile singularity --datadir data 	--input /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/sequences/traf_2000_sequences.fasta --applications Pfam,AntiFam,mobidblite --formats tsv 	--outprefix traf_2000_is6
mv traf_2000_is6.tsv /homes/ako/nfs_research/projects/foldunfold2/part3_gpc_testing/results/traf_2000_is6.tsv
