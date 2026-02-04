import os
import subprocess

OUTDIR = "part2_gpc_training/1_sequence_selection/sequences/synth_afms"
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

os.makedirs(f"{OUTDIR}/tmp_data/diamond", exist_ok=True)
os.makedirs(f"{OUTDIR}/tmp_data/diamond/logs", exist_ok=True)
os.makedirs(f"{OUTDIR}/tmp_data/is6", exist_ok=True)
os.makedirs(f"{OUTDIR}/tmp_data/is6/logs", exist_ok=True)
os.makedirs(f"{OUTDIR}/tmp_data/logs", exist_ok=True)

for org in organisms.keys():
    fasta_file = organisms[org]["fasta"]

    # filter proteome dna (translate 1 frame), then run iupred,

    subprocess.run(
        [
            "python",
            f"{SCRIPTDIR}/generate_afms.py",
            f"{REF_GENOMES}/{fasta_file}",
            org,
            "--min_length",
            "30",
        ]
    )

    subprocess.run(
        [
            "python",
            f"{SCRIPTDIR}/make_queries.py",
            f"{OUTDIR}/tmp_fasta/{org}_all_fake_seqs.fasta",
            "--resultsdir",
            f"{OUTDIR}/tmp_data/diamond/temp_queries/{org}",
        ]
    )

    with open(f"{OUTDIR}/tmp_data/logs/{org}_make_jobs.log", "w") as f:
        subprocess.run(
            [
                "python",
                f"{SCRIPTDIR}/prepare_slurm.py",
                f"{OUTDIR}/tmp_fasta/{org}_all_fake_seqs.fasta",
                org,
                "--temp_fasta_dir",
                f"{OUTDIR}/tmp_data/diamond/temp_queries/{org}",
            ],
            stdout=f,
        )

    subprocess.run(["sbatch", f"{OUTDIR}/tmp_data/diamond/{org}_run_dia.sh"])
    subprocess.run(["sbatch", f"{OUTDIR}/tmp_data/is6/{org}_run_is6.sh"])