#!/bin/bash
#SBATCH --job-name=run_checks_spaf
#SBATCH --output=part3_gpc_testing/results/tmp/run_checks_spaf.log
#SBATCH --time=2:00:00
#SBATCH --mem=500MB

python part3_gpc_testing/2_gpc_validation/job_scripts/prepare_scripts.py \\
part3_gpc_testing/1_gpc_run/results/sequences/spaf_sequences.fasta \\
spaf part3_gpc_testing/2_gpc_validation/results/spaf --tools shadow spurio

sbatch part3_gpc_testing/results/tmp/spaf_run_shadow.sh
sbatch part3_gpc_testing/results/tmp/spaf_run_spurio.sh
