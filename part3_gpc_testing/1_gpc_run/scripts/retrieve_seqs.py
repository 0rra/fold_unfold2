import argparse

import pandas as pd
from Bio import SeqIO

swissprot_sequence_path = "data/afdb_bact_only/2025_03_afdb_sprot_bact.fasta"
trembl_sequence_path = "data/afdb_bact_only/2025_03_afdb_trem_bact.fasta"
outdir = "part3_gpc_testing/1_gpc_run/results/sequences"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_pred", help="Input prediction csv e.g. spaf_predictions.csv"
    )
    parser.add_argument(
        "output_prefix", help="prefix for output fasta file, e.g. spaf_"
    )
    parser.add_argument(
        "-sample",
        "--sample",
        default="all",
        help="Number of sequences to randomly sample, or 'all' for all sequences (default: all)",
    )
    args = parser.parse_args()

    pred_csv = args.input_pred
    Out_Prefix = args.output_prefix
    sample_arg = args.sample

    source_file = ""
    if Out_Prefix.lower().startswith("sp"):
        source_file = swissprot_sequence_path
    elif Out_Prefix.lower().startswith("tr"):
        source_file = trembl_sequence_path
    else:
        print("please give a valid out prefix starting with either sp or tr")
        exit()

    ids = set()

    for batch in pd.read_csv(pred_csv, chunksize=10000, iterator=True):
        batch_ids = batch["uniprot_id"].tolist()
        ids.update(batch_ids)

    # sample sequences
    if sample_arg.lower() != "all":
        try:
            sample_size = int(sample_arg)
            if sample_size < len(ids):
                ids = set(
                    pd.Series(list(ids)).sample(n=sample_size, random_state=42).tolist()
                )
                print(f"Randomly sampled {sample_size} sequences from {len(ids)} total")
            else:
                print(
                    f"Sample size ({sample_size}) >= total sequences ({len(ids)}), using all sequences"
                )
        except ValueError:
            print(f"Invalid sample argument: {sample_arg}. Use 'all' or an integer.")
            exit()

    # use these ids to filter source_file fasta to make out_file fasta
    if sample_arg.lower() == "all":
        out_file = f"{outdir}/{Out_Prefix}_sequences.fasta"
    else:
        out_file = f"{outdir}/{Out_Prefix}_{sample_arg}_sequences.fasta"

    # filter fasta based on collected ids
    with open(out_file, "w") as out_handle:
        for record in SeqIO.parse(source_file, "fasta"):
            # checks for seq id
            if record.id in ids:
                SeqIO.write(record, out_handle, "fasta")

    print(f"Filtered {len(ids)} sequences to {out_file}")


if __name__ == "__main__":
    main()
