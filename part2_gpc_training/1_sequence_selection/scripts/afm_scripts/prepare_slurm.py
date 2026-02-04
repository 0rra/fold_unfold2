import argparse
from pathlib import Path

OUTDIR = "part2_gpc_training/1_sequence_selection/sequences/synth_afms/tmp_data"

run_is6_template = [
    "#!/bin/bash",
    "#SBATCH --job-name={prefix}_is6",
    "#SBATCH --output={OUTDIR}/is6/logs/{prefix}_is6.out",
    "#SBATCH --error={OUTDIR}/is6/logs/{prefix}_is6.err",
    "#SBATCH --time=2:00:00",
    "#SBATCH --mem=20GB",
    "",
    "module load nextflow/25.04.6"
    "nextflow run ebi-pf-team/interproscan6 -r 6.0.0 \\",
    "    -profile singularity --datadir data \\",
    "    --input {fasta} --applications Pfam,AntiFam,MobiDB-lite --formats tsv \\",
    "    --outprefix {prefix}_is6 --outdir {OUTDIR}/is6/is6_res",
]


run_dia_template = [
    "#!/bin/bash",
    "#SBATCH --job-name={prefix}_run_dia",
    "#SBATCH --output={OUTDIR}/diamond/logs/{prefix}_dia_%A_%a.out",
    "#SBATCH --error={OUTDIR}/diamond/logs/{prefix}_dia_%A_%a.err",
    "#SBATCH --time=42:00:00",
    "#SBATCH --cpus-per-task=4",
    "#SBATCH --mem=20GB",
    "#SBATCH --array=1-{jobs}%{max_jobs}",
    "",
    "source /homes/$USER/.bashrc",
    'JOB_FILE="{OUTDIR}/diamond/{prefix}_diamond_jobs.txt"',
    'DIAMOND_DB="data/ref_prot/Bacteria_ref_proteomes.dmnd"',
    'RESULTS_DIR="{OUTDIR}/diamond/diamond_res/{prefix}"',
    "THREADS=$SLURM_CPUS_PER_TASK",
    "# Load environment",
    "# Create directories",
    'mkdir -p "$RESULTS_DIR"',
    'mkdir -p "{OUTDIR}/diamond/logs"',
    "# Get the line corresponding to this array task",
    'JOB_LINE=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" "$JOB_FILE")',
    "# Check if job line exists",
    'if [ -z "$JOB_LINE" ]; then',
    '    echo "ERROR: No job found for array task $SLURM_ARRAY_TASK_ID"',
    "    exit 1",
    "fi",
    "# Convert job line to array of files",
    'read -ra FASTA_FILES <<< "$JOB_LINE"',
    "NUM_FILES=${{#FASTA_FILES[@]}}",
    'echo "=================================================="',
    'echo "Array Task ID: $SLURM_ARRAY_TASK_ID"',
    'echo "Processing $NUM_FILES files"',
    'echo "Start time: $(date)"',
    'echo "=================================================="',
    'echo ""',
    "# Counter for tracking progress",
    "PROCESSED=0",
    "SKIPPED=0",
    "FAILED=0",
    "# Process each file assigned to this job",
    'for FASTA_FILE in "${{FASTA_FILES[@]}}"; do',
    "    # Skip empty entries",
    '    if [ -z "$FASTA_FILE" ]; then',
    "        continue",
    "    fi",
    "    ",
    "    # Check if file exists",
    '    if [ ! -f "$FASTA_FILE" ]; then',
    '        echo "ERROR: File not found: $FASTA_FILE"',
    "        FAILED=$((FAILED + 1))",
    "        continue",
    "    fi",
    "    ",
    "    # Get basename for output",
    '    BASENAME=$(basename "$FASTA_FILE" .fasta)',
    '    OUTPUT_FILE="$RESULTS_DIR/${{BASENAME}}.m8"',
    "    ",
    '    echo "----------------------------------------"',
    '    echo "Processing: $BASENAME"',
    '    echo "Start: $(date)"',
    "    ",
    "    # Skip if output exists",
    '    if [ -f "$OUTPUT_FILE" ]; then',
    '        echo "Output already exists: $OUTPUT_FILE"',
    '        echo "Skipping..."',
    "        SKIPPED=$((SKIPPED + 1))",
    "        continue",
    "    fi",
    "    ",
    "    # Run DIAMOND",
    "    diamond blastp \\",
    '        -q "$FASTA_FILE" \\',
    '        -d "$DIAMOND_DB" \\',
    '        -o "$OUTPUT_FILE" \\',
    "        -f 6 sseqid qseqid qlen pident length qstart qend sstart send evalue stitle qtitle \\",
    "        -p $THREADS",
    "    ",
    "    EXIT_CODE=$?",
    "    ",
    "    if [ $EXIT_CODE -eq 0 ]; then",
    '        echo "Completed: $(date)"',
    "        PROCESSED=$((PROCESSED + 1))",
    "    else",
    '        echo "FAILED with exit code $EXIT_CODE: $(date)"',
    "        FAILED=$((FAILED + 1))",
    "    fi",
    '    echo ""',
    "done",
    'echo "=================================================="',
    'echo "Job Summary for Array Task $SLURM_ARRAY_TASK_ID"',
    'echo "Files assigned: $NUM_FILES"',
    'echo "Processed successfully: $PROCESSED"',
    'echo "Skipped (already done): $SKIPPED"',
    'echo "Failed: $FAILED"',
    'echo "End time: $(date)"',
    'echo "=================================================="',
    "# Exit with error if any files failed",
    "if [ $FAILED -gt 0 ]; then",
    "    exit 1",
    "fi",
]


def generate_job_file(fasta_dir, job_file, files_per_job=100):
    # Convert to Path objects
    fasta_dir = Path(fasta_dir)
    job_file = Path(job_file)

    job_file.parent.mkdir(parents=True, exist_ok=True)

    # Collect all fasta files
    fasta_files = sorted(fasta_dir.glob("*.fasta"))
    total_files = len(fasta_files)

    print(f"Total FASTA files found: {total_files}")

    if total_files == 0:
        print(f"Error: No FASTA files found in {fasta_dir}")
        return

    # Write batches to job file
    with open(job_file, "w") as f:
        for i in range(0, total_files, files_per_job):
            # Get batch of files
            batch = fasta_files[i : i + files_per_job]
            # Write space-separated file paths
            job_line = " ".join(str(fasta) for fasta in batch)
            f.write(job_line + "\n")

    # Calculate total jobs
    total_jobs = (total_files + files_per_job - 1) // files_per_job

    print(
        f"Generated {job_file} with {total_jobs} jobs ({files_per_job} files per job)"
    )
    return int(total_jobs)


def write_file(template, input_dict, file_name):
    # Ensure parent directory exists
    Path(file_name).parent.mkdir(parents=True, exist_ok=True)

    prefix = input_dict["prefix"]
    fasta = input_dict["fasta"]
    jobs = input_dict["jobs"]
    max_jobs = input_dict["max_jobs"]
    outdir = input_dict["OUTDIR"]

    with open(file_name, "w") as fh:
        for line in template:
            fh.write(
                line.format(
                    prefix=prefix,
                    fasta=fasta,
                    jobs=jobs,
                    max_jobs=max_jobs,
                    OUTDIR=outdir,
                )
                + "\n"
            )
    print("Created", file_name)


def main():
    parser = argparse.ArgumentParser(
        description="Writes slurm scripts for running diamond and is6"
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument(
        "--temp_fasta_dir",
        default="temp_queries",
        help="Directory path for diamond temp query fastas",
    )
    parser.add_argument("output_prefix", help="prefix for output files")
    args = parser.parse_args()

    FastaFile = args.input_fasta
    OutPrefix = args.output_prefix
    FastaDir = args.temp_fasta_dir

    JobFile = f"{OUTDIR}/diamond/{OutPrefix}_diamond_jobs.txt"

    total_jobs = generate_job_file(FastaDir, JobFile, 100)
    print(total_jobs)

    if total_jobs is None:
        print("No jobs generated. Exiting.")
        return

    if total_jobs > 10:
        max_jobs = 10
    else:
        max_jobs = total_jobs

    fill_in = {
        "prefix": OutPrefix,
        "fasta": FastaFile,
        "jobs": total_jobs,
        "max_jobs": max_jobs,
        "OUTDIR": OUTDIR,
    }

    templates = [run_is6_template, run_dia_template]
    template_names = ["is6/run_is6.sh", "diamond/run_dia.sh"]

    for template, name in zip(templates, template_names):
        dirname, filename = name.split("/")
        # scripts are written to OUTDIR/subdir/prefix_filename
        file_name = f"{OUTDIR}/{dirname}/{OutPrefix}_{filename}"
        print("Writing slurm script", file_name)
        write_file(template, fill_in, file_name)


if __name__ == "__main__":
    main()
