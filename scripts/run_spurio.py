import argparse
import os
import subprocess

from Bio import SeqIO

# to run update with path to user's spurio installation
SPURIO_PATH = "/path/to/spurio/spurio.py"


def main():
    parser = argparse.ArgumentParser(
        description="Run spurio on each sequence in a FASTA file individually."
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    args = parser.parse_args()

    FileIn = os.path.abspath(args.input_fasta)
    basename = os.path.splitext(os.path.basename(FileIn))[0]

    log_filename = os.path.abspath("{}_spurio_log.txt".format(basename))
    summary_filename = os.path.abspath("{}_summary.txt".format(basename))

    spurio_dir = os.path.dirname(SPURIO_PATH)
    print(spurio_dir)

    with open(log_filename, "w") as log_file:
        with open(FileIn) as handle:
            seq_records = SeqIO.parse(handle, "fasta")
            for i, seq_record in enumerate(seq_records):
                temp_filename = os.path.abspath(
                    "{}_temp_seq_{}.fasta".format(basename, i)
                )

                with open(temp_filename, "w") as temp_out:
                    SeqIO.write(seq_record, temp_out, "fasta")

                search_cmd = [
                    "python",
                    SPURIO_PATH,
                    "-q",
                    temp_filename,
                    "-o",
                    summary_filename,
                ]
                result = subprocess.run(
                    search_cmd,
                    cwd=spurio_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                )

                # Write everything to log
                log_file.write("Sequence {}: {}\n".format(i + 1, seq_record.id))
                log_file.write(result.stdout)
                log_file.write(result.stderr)
                log_file.write("\n---\n")

                os.remove(temp_filename)

    print("Done. Log: {}".format(log_filename))
    print("Summary: {}".format(summary_filename))


if __name__ == "__main__":
    main()
