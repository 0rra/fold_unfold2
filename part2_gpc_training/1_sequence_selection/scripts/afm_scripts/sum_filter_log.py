import glob

# summarise filter logs
start_seqs = []
blast_removed = []
pfam_removed = []
cluster_removed = []
total_removed = []
kept_seq = []

# get list of org logs, files which match filter*_UP*log
org_logs = glob.glob("filter*_UP*log")

for log in org_logs:
    with open(log) as fh:
        for line in fh:
            if line.startswith("Initial"):
                start_seqs.append(int(line.split(": ")[1].strip()))
            elif "Blast" in line:
                blast_removed.append(int(line.split(": ")[1].strip()))
            elif "Pfam" in line:
                pfam_removed.append(int(line.split(": ")[1].strip()))
            elif "cluster deduplication" in line:
                cluster_removed.append(int(line.split(": ")[1].strip()))
            elif "seqs removed" in line:
                total_removed.append(int(line.split(": ")[1].strip()))
            elif "retained entries" in line:
                kept_seq.append(int(line.split(": ")[1].strip()))

# write filter summary file
with open("filter_summary.txt", "w") as f:
    f.write(f"Initial entries: {sum(start_seqs)}\n")
    f.write(f"Removed by Blast filter: {sum(blast_removed)}\n")
    f.write(f"Removed by Pfam filter: {sum(pfam_removed)}\n")
    f.write(f"Removed by cluster deduplication: {sum(cluster_removed)}\n")
    f.write(f"Final retained entries: {sum(kept_seq)}\n")
    f.write(f"Number seqs removed: {sum(total_removed)}\n")

print(f"Processed {len(org_logs)} log files")
print("Summary written to filter_summary.txt")
