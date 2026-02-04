import os
import time

import pandas as pd
from Bio import SeqIO

SP_FASTA = "data/swiss-prot/2025_03_uniprot_sprot.fasta"
AFDB_TAX = "data/afdb_paes/tmp_afdb_info.tsv"
AFDB_FASTA = "data/afdb_paes/afdb_sequences.fasta"
PLDDT_FILE = "data/afdb_paes/2025_03_uniprot-plddt.tsv"
PTM_FILE = "data/afdb_paes/afdb_v6_bact_paes.tsv"
OUTPUT_DIR = "data/afdb_paes"


def get_sp_proteins(fasta_path):
    return {r.description.split("|")[1] for r in SeqIO.parse(fasta_path, "fasta")}


# Load AFDB tax/id/length data and update with plddt scores and review status
def load_and_merge(tax_file, plddt_file, swissprot_ids):
    df = pd.read_csv(tax_file, sep="\t")
    df["source"] = df["uniprot_id"].apply(
        lambda x: "Swiss-Prot" if x in swissprot_ids else "TrEMBL"
    )

    plddt = pd.read_csv(plddt_file, sep="\t", names=["uniprot_id", "plddt"])
    return df.merge(plddt, on="uniprot_id", how="inner")


def filter_bacteria_and_merge_ptm(df, ptm_file):
    bact_df = df[df["is_bacteria"] == "Yes"].copy()
    ptm = pd.read_csv(ptm_file, sep="\t")
    return bact_df.merge(ptm, on="afdb_id", how="inner")


def write_filtered_fasta(fasta_path, uniprot_ids, output_path):
    with open(output_path, "w") as f:
        for record in SeqIO.parse(fasta_path, "fasta"):
            uniprot_id = record.description.split("UA=")[1].split()[0]
            if uniprot_id in uniprot_ids:
                f.write(f">{uniprot_id}\n{record.seq}\n")


def main():
    start_time = time.time()

    afdb_info = OUTPUT_DIR / "2025_03_afdb_info.tsv"
    if not os.path.exists(afdb_info):
        swissprot_ids = get_sp_proteins(SP_FASTA)
        print(f"Swiss-Prot IDs extracted: {time.time() - start_time:.2f}s")

        afdb_df = load_and_merge(AFDB_TAX, PLDDT_FILE, swissprot_ids)
        afdb_df.to_csv(afdb_info, sep="\t", index=False)
        print(f"Created afdb info: {time.time() - start_time:.2f}s")

    else:
        afdb_df = pd.read_csv(afdb_info, sep="\t")
        print(f"Loaded afdb info: {time.time() - start_time:.2f}s")

    bact_info = OUTPUT_DIR / "2025_03_afdb_bact_info.tsv"
    if not os.path.exists(bact_info):
        bact_df = filter_bacteria_and_merge_ptm(afdb_df, PTM_FILE)
        bact_df.to_csv(bact_info, sep="\t", index=False)
        print(f"Created bacteria-only afdb info: {time.time() - start_time:.2f}s")

    else:
        bact_df = pd.read_csv(bact_info, sep="\t")
        print(f"Loaded bacteria-only afdb info: {time.time() - start_time:.2f}s")

    # split by review status and write sep fasta files
    sp_info = OUTPUT_DIR / "2025_03_afdb_sprot_bact_info.tsv"
    if not os.path.exists(sp_info):
        for source, label in [("Swiss-Prot", "sprot"), ("TrEMBL", "trembl")]:
            subset = bact_df[bact_df["source"] == source]
            subset.to_csv(
                OUTPUT_DIR / f"2025_03_afdb_{label}_bact_info.tsv",
                sep="\t",
                index=False,
            )
            write_filtered_fasta(
                AFDB_FASTA,
                set(subset["uniprot_id"]),
                OUTPUT_DIR / f"2025_03_afdb_{label}_bact.fasta",
            )

    print(f"\nTotal time: {time.time() - start_time:.2f}s")


if __name__ == "__main__":
    main()
