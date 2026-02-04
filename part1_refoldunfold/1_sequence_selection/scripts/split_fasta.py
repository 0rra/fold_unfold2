import sys

from Bio import SeqIO

FileIn = sys.argv[1]
des_lengths = sys.argv[2]

des_lengths = int(des_lengths)

seq_lengths = {}
with open(FileIn) as handle:
    seq_records = SeqIO.parse(handle, "fasta")
    for seq_record in seq_records:
        length = len(seq_record.seq)
        if length == des_lengths:
            print(f">{seq_record.id}\n{seq_record.seq}")
