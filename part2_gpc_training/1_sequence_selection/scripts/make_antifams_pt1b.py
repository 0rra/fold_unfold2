import os
import subprocess

OUTDIR = "part2_gpc_training/1_sequence_selection/sequences/synth_afms"

# dict of orgs
REF_GENOMES = "data/ref_prot/gen_af_proteomes"


def collect_orgs(file_path):
    organisms = {}

    for fname in os.listdir(file_path):
        if not fname.endswith(".fasta"):
            continue

        proteome_id = fname.split("_", 1)[0]

        organisms[proteome_id] = {"fasta": fname}

    return organisms


organisms = collect_orgs(REF_GENOMES)

for org in organisms.keys():
    # join all
    with open(f"{OUTDIR}/tmp_fasta/all_fake_seqs.fasta", "a") as out_f:
        with open(f"{OUTDIR}/tmp_fasta/{org}_all_fake_seqs.fasta", "r") as in_f:
            out_f.write(in_f.read())

os.makedirs(f"{OUTDIR}/tmp_data/mmseqs", exist_ok=True)
subprocess.run(
    [
        "mmseqs",
        "easy-cluster",
        f"{OUTDIR}/tmp_fasta/all_fake_seqs.fasta",
        f"{OUTDIR}/tmp_data/mmseqs/all_fake_seqs_30",
        f"{OUTDIR}/tmp_data/mmseqs/tmp",
        "--min-seq-id",
        "0.3",
    ]
)
subprocess.run(["rm", "-rf", f"{OUTDIR}/tmp_data/mmseqs/tmp"])
