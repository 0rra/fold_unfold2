import os

import pandas as pd

outdir = "part3_gpc_testing/1_gpc_run/results"

# Read file in chunks
chunks = pd.read_csv(f"{outdir}/afdb_predictions.tsv", sep="\t", chunksize=100000)
afdb = pd.concat(chunks, ignore_index=True)
print("read file")

# split fasta
spaf_file = os.path.join(outdir, "spaf_predictions.csv")
if os.path.exists(spaf_file):
    spaf_df = pd.read_csv(spaf_file, index_col=False, sep=",")
else:
    spaf_df = afdb[
        (afdb["source"] == "Swiss-Prot") & (afdb["prediction_label"] == "AntiFam")
    ]
    print(spaf_df.head())
    spaf_df.to_csv(spaf_file, index=False, sep=",")

spsp_file = os.path.join(outdir, "spsp_predictions.csv")
if os.path.exists(spsp_file):
    spsp_df = pd.read_csv(spsp_file, index_col=False, sep=",")
else:
    spsp_df = afdb[
        (afdb["source"] == "Swiss-Prot") & (afdb["prediction_label"] == "SwissProt")
    ]
    print(spsp_df.head())
    spsp_df.to_csv(spsp_file, index=False, sep=",")

traf_file = os.path.join(outdir, "traf_predictions.csv")
if os.path.exists(traf_file):
    traf_df = pd.read_csv(traf_file, index_col=False, sep=",")
else:
    traf_df = afdb[
        (afdb["source"] == "TrEMBL") & (afdb["prediction_label"] == "AntiFam")
    ]
    print(traf_df.head())
    traf_df.to_csv(traf_file, index=False, sep=",")

trsp_file = os.path.join(outdir, "trsp_predictions.csv")
if os.path.exists(trsp_file):
    trsp_df = pd.read_csv(trsp_file, index_col=False, sep=",")
else:
    trsp_df = afdb[
        (afdb["source"] == "TrEMBL") & (afdb["prediction_label"] == "SwissProt")
    ]
    print(trsp_df.head())
    trsp_df.to_csv(trsp_file, index=False, sep=",")
