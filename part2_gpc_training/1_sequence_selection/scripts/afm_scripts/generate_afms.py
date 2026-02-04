import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

OUTDIR = "part2_gpc_training/1_sequence_selection/sequences/synth_afms"


# reads in proteome dna fasta
def parse_proteome(proteome_dna_path, maxlen, minlen):
    proteome_fasta = SeqIO.to_dict(SeqIO.parse(proteome_dna_path, "fasta"))

    generated_seqs = {}
    seq_metadata = []

    frame_names = ["+2", "+3", "-1", "-2", "-3"]  # reading frames

    for gene in proteome_fasta:
        gene_seq = proteome_fasta[gene]
        try:
            # generate for 5 frames
            all_frames = [
                gene_seq.seq[1:].translate(),  # +2
                gene_seq.seq[2:].translate(),  # +3
                gene_seq.seq.reverse_complement().translate(),  # -1
                gene_seq.seq[:-1].reverse_complement().translate(),  # -2
                gene_seq.seq[:-2].reverse_complement().translate(),  # -3
            ]

            # store longest sequence per frame
            for frame_idx, entry in enumerate(all_frames):
                start_positions = [pos for pos, char in enumerate(entry) if char == "M"]

                # track longest protein in frame
                longest_seq = None
                longest_metadata = None
                max_length = minlen

                for start in start_positions:
                    # find first stop after start
                    first_stop = None
                    for pos in range(start + 1, len(entry)):
                        if entry[pos] == "*":
                            first_stop = pos
                            break

                    # option to use end of sequence if no stop
                    if first_stop is None:
                        first_stop = len(entry)
                        has_stop = False
                    else:
                        has_stop = True

                    # check if this is the longest sequence
                    seq_length = first_stop - start
                    if seq_length > max_length and seq_length < maxlen:
                        if "|" in gene:
                            gene_id = gene.split("|")[1]
                        else:
                            gene_id = gene
                        max_length = seq_length
                        longest_seq = entry[start:first_stop]
                        longest_metadata = {
                            "gene": gene,
                            "fake_protein_id": f"{gene_id}_{start}_{frame_names[frame_idx].strip('+')}",
                            "frame": frame_names[frame_idx],
                            "stop_codon_presence": has_stop,
                            "start_aa": start,
                            "stop_aa": first_stop,
                            "sequence_length": seq_length,
                        }

                # add longest protein from this frame if one was found
                if longest_seq is not None:
                    fake_protein_id = longest_metadata["fake_protein_id"]
                    generated_seqs[fake_protein_id] = longest_seq
                    seq_metadata.append(longest_metadata)

        except Exception as e:
            print(gene)
            print(e)

    return generated_seqs, seq_metadata


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

    os.makedirs(
        f"{OUTDIR}/plots",
        exist_ok=True,
    )
    output_file = f"{OUTDIR}/plots/{output_prefix}_all_length_histogram.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Generates AntiFam-like sequences from a given proteome dna fasta"
    )
    parser.add_argument(
        "input_fasta", help="Input FASTA file, must be proteome_dna.fasta format"
    )
    parser.add_argument("output_prefix", help="prefix for output file")
    parser.add_argument(
        "--min_length", help="min sequence length, default is 20", default=20
    )
    parser.add_argument(
        "--max_length", help="max sequence length, default is 1000", default=1000
    )
    args = parser.parse_args()

    FastaFile = args.input_fasta
    Out_Prefix = args.output_prefix
    MIN_LEN = int(args.min_length)
    MAX_LEN = int(args.max_length)
    print("Starting to parse proteome")
    fake_seqs, fake_seq_data = parse_proteome(FastaFile, MAX_LEN, MIN_LEN)

    fake_seq_df = pd.DataFrame(fake_seq_data)
    os.makedirs(
        f"{OUTDIR}/tmp_data",
        exist_ok=True,
    )
    fake_seq_df.to_csv(
        f"{OUTDIR}/tmp_data/{Out_Prefix}_all_fake_seq_data.csv",
        index=False,
    )
    print("Saved metadata to file", f"{Out_Prefix}_all_fake_seq_data.csv")

    # write fake_seqs to fasta file
    os.makedirs(
        f"{OUTDIR}/tmp_fasta",
        exist_ok=True,
    )
    Out_fasta = f"{OUTDIR}/tmp_fasta/{Out_Prefix}_all_fake_seqs.fasta"
    with open(Out_fasta, "w") as fasta_out:
        for protein_id, sequence in fake_seqs.items():
            fasta_out.write(f">{protein_id}\n")
            fasta_out.write(f"{sequence}\n")

    print(f"Saved sequences to file {Out_fasta}")

    print("Number of generated seqs : " + str(len(fake_seqs)))

    plot_sequence_lengths(fake_seq_df, Out_Prefix)


if __name__ == "__main__":
    main()
