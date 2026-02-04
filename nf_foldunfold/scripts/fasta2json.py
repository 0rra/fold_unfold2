import json
import os
import sys


def main():
    args = sys.argv[1:]
    try:
        parse_fasta(args[0], args[1], args[2])
    except IndexError:
        args.append(os.getcwd())
        parse_fasta(args[0], args[1], args[2])


def parse_fasta(fasta_file, output_dir, use_msa):
    if output_dir.endswith("/"):
        output_dir = output_dir[:-1]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print("writing to", output_dir)
    with open(fasta_file, "r") as f:
        seqs = {}
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seqid = line[1:]
                seqs[seqid] = ""
            if not line.startswith(">"):
                seqs[seqid] += line.upper()

        for seq in seqs:
            inputs = {
                "name": seq,
                "sequences": [],
                "modelSeeds": [1],
                "dialect": "alphafold3",
                "version": 1,
            }
            inputs["sequences"] = [{"protein": {"id": ["A"], "sequence": seqs[seq]}}]
            if str(use_msa) == "false":
                inputs["sequences"][0]["protein"]["unpairedMsa"] = ""
                inputs["sequences"][0]["protein"]["pairedMsa"] = ""

            with open(f"{output_dir}/{seq}.json", "w") as jf:
                jf.write(json.dumps(inputs))


if __name__ == "__main__":
    main()
