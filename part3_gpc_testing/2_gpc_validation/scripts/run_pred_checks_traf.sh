#!/bin/bash
#SBATCH --job-name=run_checks_traf
#SBATCH --output=part3_gpc_testing/results/tmp/run_checks_traf.log
#SBATCH --time=2:00:00
#SBATCH --mem=500MB

python part3_gpc_testing/2_gpc_validation/job_scripts/prepare_scripts.py \\
part3_gpc_testing/1_gpc_run/results/sequences/traf_1000_sequences.fasta \\
traf_1000 part3_gpc_testing/2_gpc_validation/results/traf_1000

sbatch part3_gpc_testing/results/tmp/traf_1000_run_shadow.sh
sbatch part3_gpc_testing/results/tmp/traf_1000_run_is6.sh
sbatch part3_gpc_testing/results/tmp/traf_1000_run_spurio.sh