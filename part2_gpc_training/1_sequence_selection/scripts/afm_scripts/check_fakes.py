import argparse
import glob
import os

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

# script to filter all fakes fasta and update metadata csv using is6 and diamond results
OUTDIR = "part2_gpc_training/1_sequence_selection/sequences/synth_afms"


def plot_sequence_lengths(seq_metadata, output_prefix):
    lengths = seq_metadata["sequence_length"].tolist()

    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=100, edgecolor="black", alpha=0.7)
    plt.xlabel("Sequence Length", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.title("Distribution of Generated Protein Sequence Lengths", fontsize=14)
    plt.grid(axis="y", alpha=0.3)

    mean_len = sum(lengths) / len(lengths)
    plt.axvline(
        mean_len,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean: {mean_len:.1f}",
    )
    plt.legend()

    output_file = f"{OUTDIR}/plots/{output_prefix}_filtered_length_histogram.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def filter_fasta(source_file, no_match_ids, out_fasta):
    with open(out_fasta, "w") as out_handle:
        for record in SeqIO.parse(source_file, "fasta"):
            if record.id in no_match_ids:
                SeqIO.write(record, out_handle, "fasta")
    return


def anno_is6(is6_path, result_df):
    is6 = pd.read_csv(is6_path, sep="\t", index_col=False, header=None)

    # make cols and split accession col
    header = [
        "UniProt ID",
        "md5",
        "Sequence length",
        "DB",
        "Family ID",
        "Desc",
        "Start",
        "End",
        "Evalue",
        "2",
        "3",
        "InterPro ID",
        "IP name",
        "4",
        "5",
    ]
    # add header
    is6.columns = header

    is6 = is6[["UniProt ID", "DB", "Family ID", "Desc", "Start", "End"]]

    ids = result_df[["fake_protein_id"]].drop_duplicates()

    # count pfam/antifam per id
    counts = (
        is6.groupby(["UniProt ID", "DB"]).size().unstack(fill_value=0).reset_index()
    )

    # merge back to keep all ids, even if counts are missing
    summary = ids.merge(
        counts, right_on="UniProt ID", left_on="fake_protein_id", how="left"
    ).fillna(0)

    # make sure both columns exist
    for col in ["Pfam", "AntiFam", "MobiDB-lite"]:
        if col not in summary.columns:
            summary[col] = 0

    # clean column names
    summary = summary.rename(
        columns={
            "Pfam": "Pfam count",
            "AntiFam": "AntiFam count",
            "MobiDB-lite": "MobiDB count",
        }
    )
    summary = summary.drop(columns=["UniProt ID"])

    result_df = result_df.merge(summary, how="left", on="fake_protein_id")

    result_df["IS6 result"] = ~(
        (result_df["Pfam count"] == 0)
        & (result_df["AntiFam count"] == 0)
        & (result_df["MobiDB count"] == 0)
    )

    return is6, result_df


def parse_m8(m8_file):
    # sseqid qseqid qlen pident length qstart qend sstart send evalue stitle qtitle
    with open(m8_file, "r") as fh:
        total_count = 0
        for line in fh:
            total_count += 1

    return total_count


def anno_dia(match_data, diamond_res):
    dia_res_files = glob.glob(f"{diamond_res}/*.m8")
    protein_matches = {}

    for file in dia_res_files:
        hit_count = 0
        id = file.replace(f"{diamond_res}/query_", "").strip(".m8")
        if os.path.getsize(file) > 0:
            hit_count = parse_m8(file)
        protein_matches[id] = hit_count
    # print(protein_matches)
    # map onto match data using fake_protein_id to make new col blast matches
    match_data["Blast_matches"] = match_data["fake_protein_id"].map(protein_matches)
    # print(match_data)
    return match_data


def check_clusters(metadata, cluster_res):
    clusters = {}
    cluster_set = set()
    with open(cluster_res) as cr:
        for line in cr:
            cluster, seq = line.strip().split("\t")
            clusters[seq] = cluster
            cluster_set.add(cluster)
    print("TOTAL CLUSTERS", len(cluster_set))
    # then map clusters to match_data
    metadata["Cluster"] = metadata["fake_protein_id"].map(clusters)
    print("UNIQUE CLUSTERS FOR ORG", len(metadata["Cluster"].unique()))
    return metadata


def filter_fakes(match_data):
    # print(match_data)
    initial_n = len(match_data)
    print(f"Initial entries: {initial_n}")

    # ---- Filter 1: Blast matches == 0 ----
    step1 = match_data[match_data["Blast_matches"] == 0]
    removed_1 = initial_n - len(step1)
    print(f"Removed by Blast filter: {removed_1}")

    # ---- Filter 2: Pfam count == 0 ----
    step2 = step1[step1["Pfam count"] == 0]
    removed_2 = len(step1) - len(step2)
    print(f"Removed by Pfam filter: {removed_2}")

    # ---- Filter 3: unique clusters ----
    step3 = step2.drop_duplicates(subset="Cluster", keep="first")
    removed_3 = len(step2) - len(step3)
    print(f"Removed by cluster deduplication: {removed_3}")

    print(f"Final retained entries: {len(step3)}")

    filtered_match_ids = step3["fake_protein_id"].tolist()

    return filtered_match_ids, step3


def main():
    parser = argparse.ArgumentParser(
        description="Checks diamond and is6 data to return fake fakes"
    )
    parser.add_argument("input_fasta", help="Path to all fakes sequences")
    parser.add_argument("output_prefix", help="Prefix for output filtered fakes")
    args = parser.parse_args()

    FastaFile = args.input_fasta
    Out_Prefix = args.output_prefix

    MetaData = f"{OUTDIR}/tmp_data/{Out_Prefix}_all_fake_seq_data.csv"
    DiaResPath = f"{OUTDIR}/tmp_data/diamond/diamond_res/{Out_Prefix}"
    IS6ResPath = f"{OUTDIR}/tmp_data/is6/is6_res/{Out_Prefix}_is6.tsv"

    metadata_df = pd.read_csv(MetaData)
    metadata_df = anno_dia(metadata_df, DiaResPath)
    _, metadata_df = anno_is6(IS6ResPath, metadata_df)

    # sort by blast hits
    metadata_df = metadata_df.sort_values(by="Blast_matches", ascending=False)

    CLUSTER_PATH = f"{OUTDIR}/tmp_data/mmseqs/all_fake_seqs_30_cluster.tsv"
    # adds sequence clusters
    metadata_df = check_clusters(metadata_df, CLUSTER_PATH)

    metadata_df["Organism"] = Out_Prefix

    # filter to get fake fakes
    filtered_match_ids, filtered_metadata = filter_fakes(metadata_df)
    # add to metadata, failed or not if in filtered_metadata
    metadata_df["Fake"] = metadata_df["fake_protein_id"].isin(filtered_match_ids)
    metadata_df.to_csv(f"{OUTDIR}/tmp_data/{Out_Prefix}_anno_all_fake_seq_data.tsv", sep="\t")

    # write fake_seqs to fasta file
    # get all fakes
    all_fake_seqs = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

    Out_fasta = f"{OUTDIR}/tmp_fasta/{Out_Prefix}_fake_fake_seqs.fasta"
    with open(Out_fasta, "w") as fasta_out:
        for protein_id, sequence in all_fake_seqs.items():
            if protein_id in filtered_match_ids:
                fasta_out.write(f">{protein_id}\n")
                fasta_out.write(f"{sequence.seq}\n")

    print(f"Saved sequences to file {Out_fasta}")

    print("Number seqs removed: ", len(metadata_df) - len(filtered_match_ids))
    print("Number of fake seqs pass filter: " + str(len(filtered_match_ids)))

    plot_sequence_lengths(filtered_metadata, Out_Prefix)


if __name__ == "__main__":
    main()