import argparse
from pathlib import Path

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", help="Multi protein fasta file")
    parser.add_argument("--resultsdir", help="optional path to results directory")
    args = parser.parse_args()

    fasta_seqs = SeqIO.to_dict(SeqIO.parse(args.fasta_file, "fasta"))
    print(f"splitting {args.fasta_file} into {len(fasta_seqs)} temp query fasta files")

    resdir = Path(args.resultsdir)
    resdir.mkdir(parents=True, exist_ok=True)

    for prot in fasta_seqs:
        prot_seq = fasta_seqs[prot]
        Out_fasta = f"{resdir}/query_{prot}.fasta"
        with open(Out_fasta, "w") as fasta_out:
            fasta_out.write(f">{prot_seq.id}\n")
            fasta_out.write(f"{prot_seq.seq}\n")

    print("Complete")


if __name__ == "__main__":
    main()
