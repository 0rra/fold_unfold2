from Bio import SeqIO
import sys

FileIn = sys.argv[1]

with open(FileIn) as handle:
    seq_records = SeqIO.parse(handle, "fasta")
    for seq_record in seq_records:
        id = seq_record.id
        name = id.split("|")[1]
        print(f">{name}")
        print(f"{seq_record.seq}")
