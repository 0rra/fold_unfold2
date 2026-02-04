import math
import os

import matplotlib.pyplot as plt
import pandas as pd

OUTDIR = "part2_gpc_training/1_sequence_selection/sequences/synth_afms"

# selection 200 from each pool of sequences
# read metadata, filtered seqs

REF_GENOMES = "data/ref_prot/gen_af_proteomes"


MIN_LENGTH = 30


def plot_sequence_lengths(seq_dict, output_prefix):
    lengths = [len(seq) for seq in seq_dict.values()]

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

    # Save plot
    output_file = f"{OUTDIR}/plots/{output_prefix}_length_histogram.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Plotted {len(lengths)} sequences")
    print(f"Mean length: {mean_len:.1f}")
    print(f"Saved plot to {output_file}")


def sample_seqs(out_fasta, org):
    fasta_dict = {}
    with open(out_fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                id = line.strip(">").strip()
                fasta_dict[id] = ""
            else:
                fasta_dict[id] += line.strip()
    # filter very short!
    print("NUMBER TOTAL SEQS", len(fasta_dict), org)
    longer_fasta_dict = {
        key: fasta_dict[key]
        for key in fasta_dict.keys()
        if len(fasta_dict[key]) > MIN_LENGTH
    }
    # sample across bins
    # make a dataframe, sequence id, sequence length
    print("NUMBER PASS MIN LENGTH", len(longer_fasta_dict), org)
    id_lengths = {id: len(longer_fasta_dict[id]) for id in longer_fasta_dict.keys()}
    # then make this a dataframe
    lengths_df = pd.DataFrame(list(id_lengths.items()), columns=["id", "length"])

    sample_size = 150
    bins = [25, 50, 75, 100, 120, 140, 160, 180, 200, 400, 600, 800, 1000]
    n_per_bin = math.ceil(sample_size / (len(bins) - 1))
    lengths_df["length_bin"] = pd.cut(lengths_df["length"], bins=bins)

    # now sample equally from each bin
    sampled_df = (
        lengths_df.groupby("length_bin", observed=True)
        .apply(
            lambda x: x.sample(min(len(x), n_per_bin), random_state=42),
            include_groups=False,
        )
        .reset_index(drop=True)
    )

    sampled_ids = sampled_df["id"].tolist()
    print("SAMPLED", len(sampled_ids), org)

    sampled_dict = {key: longer_fasta_dict[key] for key in sampled_ids}

    return sampled_dict


def collect_orgs(file_path):
    organisms = {}

    for fname in os.listdir(file_path):
        if not fname.endswith(".fasta"):
            continue

        proteome_id = fname.split("_", 1)[0]

        organisms[proteome_id] = {"fasta": fname}

    return organisms


organisms = collect_orgs(REF_GENOMES)
combined_sampled = {}
metadata_dfs = []
for org in organisms:
    fasta_file = f"{OUTDIR}/tmp_fasta/{org}_fake_fake_seqs.fasta"
    sampled = sample_seqs(fasta_file, org)
    combined_sampled.update(sampled)
    metadata_file = f"{OUTDIR}/tmp_data/{org}_anno_all_fake_seq_data.tsv"
    metadata_df = pd.read_csv(metadata_file, sep="\t", index_col=0)
    metadata_dfs.append(metadata_df)

combined_metadata = pd.concat(metadata_dfs)

print(f"{len(combined_sampled)} sequences selected! ")

# then filter to make metadata file of sampled sequences
combined_metadata.to_csv(f"{OUTDIR}/tmp_data/combined_anno_fake_seq_data.tsv", sep="\t")
print(combined_metadata.head())
filtered_metadata = combined_metadata[
    combined_metadata["fake_protein_id"].isin(combined_sampled.keys())
]

# check no cluster matches:
original_length = len(filtered_metadata)
filtered_metadata = filtered_metadata.drop_duplicates(subset="Cluster", keep="first")
if len(filtered_metadata) != original_length:
    print(
        f"Warning: Removed {original_length - len(filtered_metadata)} non-unique clusters"
    )
    print(f"Final count: {len(filtered_metadata)} unique sequences")

filtered_metadata.to_csv(f"{OUTDIR}/tmp_data/sampled_anno_fake_seq_data.tsv", sep="\t")

# Update combined_sampled to only include sequences in filtered_metadata
combined_sampled = {
    key: combined_sampled[key] for key in filtered_metadata["fake_protein_id"]
}

print(
    "WARNING PFAM MATCHES: ",
    len(filtered_metadata[filtered_metadata["Pfam count"] > 0]),
)
print(
    "WARNING DIAMOND MATCHES: ",
    len(filtered_metadata[filtered_metadata["Blast_matches"] > 0]),
)

plot_sequence_lengths(combined_sampled, "sampled_fakes")

# Write to file
with open(f"{OUTDIR}/sampled_fakes_orgs.fasta", "w") as fh:
    for seq in combined_sampled.keys():
        fh.write(f">{seq}\n{combined_sampled[seq]}\n")
