#!/bin/bash
out_dir="part3_gpc_testing/1_gpc_run/results"
script_dir="part3_gpc_testing/1_gpc_run/scripts"

python ${script_dir}/split_preds.py

python ${script_dir}/retrieve_seqs.py ${out_dir}/spaf_predictions.csv spaf -sample all

python ${script_dir}/retrieve_seqs.py ${out_dir}/traf_predictions.csv traf -sample 1000