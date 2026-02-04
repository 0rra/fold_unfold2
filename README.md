# FoldUnfold2

## Overview

Code to reproduce analysis in "Folding the unfoldable 2: using AlphaFold and ESMFold to explore spurious proteins".

## Project Structure
```
foldunfold2/
├── data/                    # Data preparation
├── part1_refoldunfold/      # Investigation of confidence scores across structure prediction methods and sequence types
├── part2_gpc_training/      # GPC model training and evaluation
└── part3_afdb_gpc/          # Application of model to AlphaFold Database proteins
```

---

## Data Sources

- **AntiFam**: Seed sequences from [AntiFam](https://github.com/ebi-pf-team/antifam)
- **SwissProt**: Sequences from [UniProt](https://www.uniprot.org/) (release 2025_03)
- **AFDB**: Structures and confidence scores from [AlphaFold Database](https://alphafold.ebi.ac.uk/) (v6)
- **Proteomes**: Reference proteomes from [UniProt](https://www.uniprot.org/) (release 2025_03)

---

## Software dependencies

- **MMseqs2**: sequence clustering ([GitHub](https://github.com/soedinglab/MMseqs2))
- **IUPred2a**: disorder prediction ([website](https://iupred2a.elte.hu))
- **InterProScan6**: protein family and functional annotation ([GitHub](https://github.com/ebi-pf-team/interproscan6))
- **Diamond**: protein sequence search ([GitHub](https://github.com/bbuchfink/diamond))
- **Spurio**: spurious ORF detection ([Bitbucket](https://bitbucket.org/bateman-group/spurio))
---

## Environment Set Up

```bash
conda env create -f environment.yml
conda activate foldunfold2
```

---

## Structure Prediction Pipeline

All structure predictions use the Nextflow pipeline in `nf_foldunfold/` (nextflow version 25.04.6):
```bash
nextflow run nf_foldunfold/main.nf \
    --predictmethods colabfold,esmfold,alphafold3 \
    --input inputs/sequences.fasta \
    --outdir output/structure_prediction_dir

# Parse results
python scripts/get_res.py output/structure_prediction_dir output/parsed_results --prefix res_ 
```

---

## Part 1: FoldUnfold Recreation

Comparing structure prediction confidence between AntiFam (spurious), SwissProt (real), and random proteins.

### Sequence Selection
```bash
python part1_refoldunfold/1_sequence_selection/scripts/seed2seq.py
python part1_refoldunfold/1_sequence_selection/scripts/make_random_seqs.py \
    10,16,20,30,40,50,60,70,80,90,100,120,140,160,180,200 5 swissprot > part1_refoldunfold/1_sequence_selection/random_seqs/weighted_rdm.fasta

./part1_refoldunfold/1_sequence_selection/scripts/pick_sp.sh
```

### Structure prediction results
Parsed confidence scores for AntiFam, random and Swiss-Prot structure predictions
```
part1_refoldunfold/2_structure_predictions/parsed_results/
```

### Results
```bash
part1_refoldunfold/3_results_plotting/plot_figures.ipynb
```

---

## Part 2: GPC Training

Training Gaussian Process Classifier using structure prediction confidence scores (length, pTM, mean pLDDT).

### Swiss-Prot sequence Selection
```bash
python part2_gpc_training/1_sequence_selection/scripts/sample_swissprots.py
```

### Pseudo-AntiFam Generation

The synthetic antifam generation scripts were written to work with SLURM. 

```bash
python part2_gpc_training/1_sequence_selection/scripts/make_antifams_pt1.py
python part2_gpc_training/1_sequence_selection/scripts/make_antifams_pt1b.py
python part2_gpc_training/1_sequence_selection/scripts/make_antifams_pt2.py
python part2_gpc_training/1_sequence_selection/scripts/make_antifams_pt3.py
```

### Structure prediction results
Parsed confidence scores for Swiss-Prot structure predictions used for GPC testing and training. 
```
part2_gpc_training/model_data/parsed_results/
```

### Training & Evaluation
```bash
part2_gpc_training/3_train_test/gpc_dev.ipynb
```

---

## Part 3: AFDB Application

Applying GPC to bacterial proteins in AlphaFold DB.

### Prediction & Validation
```bash
# model predictions using alphafolddb features
python part3_afdb_gpc/1_gpc_run/scripts/afdb_search.py
# investigate proteins GPC predicts as spurious
./part3_afdb_gpc/1_gpc_run/scripts/process_preds.sh
python part3_gpc_testing/2_gpc_validation/scripts/prepare_scripts.py
./part3_afdb_gpc/2_gpc_validation/scripts/run_pred_checks_spaf.sh   # Spurious predicted Swiss-Prot proteins
./part3_afdb_gpc/2_gpc_validation/scripts/run_pred_checks_traf.sh   # Spurious predicted TrEMBL proteins
```

### Results
```bash
part3_gpc_testing/3_results_plotting/sum_preds.ipynb
```

