import csv

from Bio import SeqIO

# Load bacterial taxids
tax_file = "data/taxonomy/taxonomy_ancestor_2_2025_10_17.tsv"
bac_taxids = set()
with open(tax_file) as f:
    for line in f:
        if line.startswith("Taxon"):
            continue
        parts = line.split("\t")
        lineage = parts[1].split(", ")
        if len(lineage) > 1 and "Bacteria" in lineage[1]:
            bac_taxids.add(parts[0])

# Parse FASTA once, write both outputs
fasta_path = "data/afdb_paes/afdb_sequences.fasta"
with (
    open("data/afdb_paes/tmp_afdb_info.tsv", "w", newline="") as info_f,
    open("bact_afdb_ids.tsv", "w", newline="") as bact_f,
):
    info_writer = csv.writer(info_f, delimiter="\t")
    bact_writer = csv.writer(bact_f, delimiter="\t")

    info_writer.writerow(["afdb_id", "uniprot_id", "taxid", "is_bacteria", "length"])
    bact_writer.writerow(["afdb_id", "uniprot_id"])

    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description
        afdb_id = header.split(":")[1].split()[0]
        uniprot_id = header.split("UA=")[1].split()[0]
        ox = header.split("OX=")[1].split()[0]
        is_bacteria = ox in bac_taxids

        info_writer.writerow(
            [afdb_id, uniprot_id, ox, "Yes" if is_bacteria else "No", len(record.seq)]
        )
        if is_bacteria:
            bact_writer.writerow([afdb_id, uniprot_id])
