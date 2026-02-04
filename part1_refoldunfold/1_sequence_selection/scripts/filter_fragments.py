# reads fasta file and outputs new fasta file without fragments

import requests
from Bio import SeqIO
import sys
import time

FileIn = sys.argv[1]

with open(FileIn) as handle:
    seq_records = SeqIO.parse(handle, "fasta")
    for seq_record in seq_records:
        r= requests.get(f"https://rest.uniprot.org/uniprotkb/{seq_record.id}")
        info = r.json()
        try:
            if info["proteinDescription"]["flag"] == "fragment":
                pass
        except KeyError:
            print(f">{seq_record.id}\n{seq_record.seq}")
        time.sleep(0.5)
