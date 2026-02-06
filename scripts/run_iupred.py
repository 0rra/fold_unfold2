import argparse
import os
import subprocess
import sys

from Bio import SeqIO

# to run update with path to user's iupred installation
IUPRED_PATH = "/path/to/iupred2a/iupred2a.py"
iupred_dir = os.path.dirname(IUPRED_PATH)
print(iupred_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Run IUPred2A on each sequence in a FASTA file individually."
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_file", help="Output file to store IUPred2A results")
    args = parser.parse_args()

    FileIn = args.input_fasta
    OutputFile = args.output_file

    with open(FileIn) as handle, open(OutputFile, "w") as out_handle:
        seq_records = SeqIO.parse(handle, "fasta")

        for i, seq_record in enumerate(seq_records):
            temp_filename = f"{FileIn}_temp_seq_{i}.fasta"
            with open(temp_filename, "w") as temp_out:
                SeqIO.write(seq_record, temp_out, "fasta")
            # can update with path to your iupred2a dir
            search_cmd = ["python", IUPRED_PATH, temp_filename, "long"]
            result = subprocess.run(
                search_cmd,
                cwd=iupred_dir,
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                print(f"Error running iupred2a.py on {seq_record.id}", file=sys.stderr)
                print(result.stderr, file=sys.stderr)
            else:
                out_handle.write(f"# Results for {seq_record.id}\n")
                out_handle.write(result.stdout + "\n")

            os.remove(temp_filename)


if __name__ == "__main__":
    main()
