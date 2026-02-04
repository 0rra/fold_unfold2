import argparse

run_is6_template = [
    "#!/bin/bash",
    "#SBATCH --job-name={prefix}_is6",
    "#SBATCH --output={prefix}_is6.log",
    "#SBATCH --time=24:00:00",
    "#SBATCH --mem=48GB",
    "",
    "module load nextflow/25.04.6"
    "nextflow run ebi-pf-team/interproscan6 -r 6.0.0 \
	-profile singularity --datadir data \
	--input {fasta} --applications Pfam,AntiFam,mobidblite --formats tsv \
	--outprefix {prefix}_is6",
    "mv {prefix}_is6.tsv {resultsdir}/{prefix}_is6.tsv",
]

run_shadow_template = [
    "#!/bin/bash",
    "#SBATCH --job-name={prefix}_shadow",
    "#SBATCH --output={prefix}_shadow.out",
    "#SBATCH --error={prefix}_shadow.err",
    "#SBATCH --time=45:00:00",
    "#SBATCH --mem=12GB",
    "source /homes/$USER/.bashrc",
    "conda activate foldunfold",
    "cd path/to/shadow/dir"
    "python -u shadow_test.py {fasta} --resultsdir {resultsdir}/shadow_res --mode diamond --csv --trembl --bacteria",
]

run_spurio_template = [
    "#!/bin/bash",
    "#SBATCH --job-name={prefix}_spurio",
    "#SBATCH --output={prefix}_spurio.out",
    "#SBATCH --error={prefix}_spurio.err",
    "#SBATCH --time=48:00:00",
    "#SBATCH --mem=12GB",
    "",
    "source /homes/${USER}/.bashrc",
    "conda activate spurio_env2",
    "python -u scripts/run_spurio.py {fasta}",
    "mv output/summaries/query_summary.txt output/summaries/{prefix}_summary.txt",
    "cp output/summaries/{prefix}_summary.txt {resultsdir}/{prefix}_summary.txt",
]


def write_file(template, input_dict, file_name):
    prefix = input_dict["prefix"]
    fasta = input_dict["fasta"]
    resultsdir = input_dict["resultsdir"]
    with open(file_name, "w") as fh:
        for line in template:
            fh.write(
                line.format(prefix=prefix, fasta=fasta, resultsdir=resultsdir) + "\n"
            )
    print("Created", file_name)


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("input_fasta", help="Input full path FASTA file")
    parser.add_argument("output_prefix", help="prefix for output files")
    parser.add_argument("results_dir", help="directory for results")
    parser.add_argument(
        "--tools",
        nargs="+",
        choices=["shadow", "is6", "spurio"],
        default=["shadow", "is6", "spurio"],
        help="Tools to run (default: all three)",
    )
    args = parser.parse_args()

    FastaFile = args.input_fasta
    Out_Prefix = args.output_prefix
    Resdir = args.results_dir

    fill_in = {"prefix": Out_Prefix, "fasta": FastaFile, "resultsdir": Resdir}

    # Map tool names to templates
    tool_map = {
        "shadow": (run_shadow_template, "run_shadow.sh"),
        "is6": (run_is6_template, "run_is6.sh"),
        "spurio": (run_spurio_template, "run_spurio.sh"),
    }

    # Only process selected tools
    for tool in args.tools:
        template, script_name = tool_map[tool]
        file_name = Out_Prefix + "_" + script_name
        # make slurm scripts
        write_file(template, fill_in, file_name)
        print(f"Generated: {file_name}")


if __name__ == "__main__":
    main()
