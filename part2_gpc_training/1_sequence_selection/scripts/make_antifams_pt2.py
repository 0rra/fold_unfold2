import os
import subprocess

OUTDIR = "part2_gpc_training/1_sequence_selection/sequences"
SCRIPTDIR = "part2_gpc_training/1_sequence_selection/scripts/afm_scripts"
# dict of orgs
REF_GENOMES = "data/ref_prot/gen_af_proteomes"


def collect_orgs(file_path):
    organisms = {}

    for fname in os.listdir(file_path):
        if not fname.endswith(".fasta"):
            continue

        # Split on "_" and take the proteome ID
        proteome_id = fname.split("_", 1)[0]

        organisms[proteome_id] = {"fasta": fname}

    return organisms


organisms = collect_orgs(REF_GENOMES)

os.makedirs(f"{OUTDIR}/tmp_data/logs", exist_ok=True)

for org in organisms.keys():
    fasta_file = organisms[org]["fasta"]

    with open(f"{OUTDIR}/tmp_data/logs/filter_{org}.log", "w") as f:
        subprocess.run(
            [
                "python",
                f"{SCRIPTDIR}/check_fakes.py",
                f"{OUTDIR}/tmp_fasta/{org}_all_fake_seqs.fasta",
                org,
            ],
            stdout=f,
        )
