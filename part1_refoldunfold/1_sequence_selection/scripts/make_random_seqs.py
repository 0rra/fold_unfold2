import random
import sys

lengths = sys.argv[1]
seq_num = sys.argv[2]
try:
    lengths = lengths.split(",")
except ValueError:
    lengths = [lengths]
lengths = [int(x) for x in lengths]
seq_num = int(seq_num)

try:
    weight_mode = sys.argv[3]
except IndexError:
    weight_mode = "swissprot"


if weight_mode not in ["swissprot", "even"]:
    print(
        "Warning {weight_mode} not accepted, use swissprot or even\nRunning swissprot by default"
    )
    exit()

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


def protein(length, weight_mode):
    if weight_mode == "swissprot":
        return "".join(
            random.choices(
                "ARNDCQEGHILKMFPSTWYV", weights=COMPOSITION.values(), k=length
            )
        )
    else:
        return "".join(random.choices("ARNDCQEGHILKMFPSTWYV", k=length))


def make_proteins(length, weight_mode):
    for i in range(seq_num):
        print(f">{length}_{i}")
        print(protein(length, weight_mode))


for i in lengths:
    make_proteins(i, weight_mode)
