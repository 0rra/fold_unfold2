import random

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

BINS = list(range(0, 201, 20))
RANDOM_SEED = 42
ANTIFAM_FASTA = "part1_refoldunfold/1_sequence_selection/sequences/antifam_seed_seqs/seed_seqs.fasta"
SWISSPROT_TSV = "data/swiss-prot/2025_03_uniprot_sprot_anno2.tsv"


# weights obtained from https://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html
COMPOSITION = {
    "Ala": 8.25,
    "Arg": 5.53,
    "Asn": 4.06,
    "Asp": 5.45,
    "Cys": 1.37,
    "Gln": 3.93,
    "Glu": 6.75,
    "Gly": 7.07,
    "His": 2.27,
    "Ile": 5.96,
    "Leu": 9.66,
    "Lys": 5.84,
    "Met": 2.42,
    "Phe": 3.86,
    "Pro": 4.70,
    "Ser": 6.56,
    "Thr": 5.34,
    "Trp": 1.08,
    "Tyr": 2.92,
    "Val": 6.87,
}


# get AntiFam bin counts

antifam_df = pd.DataFrame(
    [{"length": len(r.seq)} for r in SeqIO.parse(ANTIFAM_FASTA, "fasta")]
)
antifam_df = antifam_df[antifam_df["length"] <= BINS[-1]]
antifam_df["length_bin_start"] = (antifam_df["length"] // 20) * 20
antifam_bin_counts = antifam_df.groupby("length_bin_start").size()

print("AntiFam counts per bin:")
for b, count in antifam_bin_counts.items():
    print(f"  {b}-{b + 20}: {count}")


# sample SwissProt seqs

df = pd.read_csv(SWISSPROT_TSV, sep="\t")

df = df[
    (df["is_fragment"] == "No")
    & (df["iupred_average"] <= 0.5)
    & (df["AntiFam count"] == 0)
]

df = df[df["length"].between(BINS[0], BINS[-1])]
df["length_bin_start"] = (df["length"] // 20) * 20

sampled_frames = []
for bin_start, group in df.groupby("length_bin_start"):
    n = antifam_bin_counts.get(bin_start, 0)
    if n == 0:
        print(f"Warning: no AntiFam seqs in bin {bin_start}-{bin_start + 20}, skipping")
        continue
    if len(group) < n:
        print(
            f"Warning: bin {bin_start}-{bin_start + 20} has only {len(group)} seqs, requested {n}"
        )
        n = len(group)
    sampled_frames.append(group.sample(n, random_state=RANDOM_SEED))

sampled = pd.concat(sampled_frames)
sampled_ids = sampled["uniprot_id"].tolist()
print(f"\nSampled {len(sampled_ids)} SwissProt sequences total")

sampled.to_csv(
    "part1_refoldunfold/1_sequence_selection/sequences/swissprot_seqs/sp_sampled.tsv",
    sep="\t",
    index=False,
)

sampled = pd.concat(sampled_frames)
sampled_ids = sampled["uniprot_id"].tolist()

# SP bin counts
print("\nSP sampled counts per bin:")
sp_sampled_bin_counts = sampled.groupby("length_bin_start").size()
for b, count in sp_sampled_bin_counts.items():
    print(f"  {b}-{b + 20}: {count}")
print(f"  Total: {len(sampled)}")

# filter and write SwissProt seqs


def filter_seqs(id_list, fasta_file, output_file):
    filtered = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        accession = record.id.split("|")[1] if "|" in record.id else record.id
        if accession in id_list:
            filtered.append(record)
    SeqIO.write(filtered, output_file, "fasta")
    return filtered


filter_seqs(
    sampled_ids,
    "data/swiss-prot/2025_03_uniprot_sprot.fasta",
    "part1_refoldunfold/1_sequence_selection/sequences/swissprot_seqs/sp_seqs.fasta",
)


# generate random sequences per bin matching AntiFam counts

random.seed(RANDOM_SEED)


def random_protein(length):
    return "".join(
        random.choices("ARNDCQEGHILKMFPSTWYV", weights=COMPOSITION.values(), k=length)
    )


random_records = []
for bin_start, count in antifam_bin_counts.items():
    bin_end = bin_start + 20
    for i in range(count):
        length = random.randint(bin_start if bin_start > 0 else 1, bin_end - 1)
        seq = random_protein(length)
        record = SeqIO.SeqRecord(
            seq=Seq(seq),
            id=f"random_{bin_start}_{bin_end}_{i}",
            description="",
        )
        random_records.append(record)

SeqIO.write(
    random_records,
    "part1_refoldunfold/1_sequence_selection/sequences/random_seqs/random_seqs2.fasta",
    "fasta",
)
print(f"Generated {len(random_records)} random sequences")

print("\nRandom sequence counts per bin:")
random_df = pd.DataFrame(
    [{"length_bin_start": int(r.id.split("_")[1])} for r in random_records]
)
random_bin_counts = random_df.groupby("length_bin_start").size()
for b, count in random_bin_counts.items():
    print(f"  {b}-{b + 20}: {count}")
print(f"  Total: {len(random_records)}")
