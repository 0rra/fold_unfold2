import os

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

DATADIR = "data/swiss-prot/"
SEQPATH = "data/swiss-prot/2025_03_uniprot_sprot_head_fixed.fasta"

SP_PATH = "part2_gpc_training/model_data/sequences"

CONFIG = {
    "swiss_prot_res": f"{DATADIR}/2025_03_uniprot_sprot_anno_afdb.tsv",
    "feats": ["pae pTM", "pLDDT mean", "length"],
    "pos_label": "AntiFam",
    "random_state": 42,
}


TEST_TRAIN_COUNTS = {
    "train": 277,
    "test50": 300,
    "sp_test98": 4900,
}


def select_seqs(ids, fasta_file, outfile):
    all_seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    with open(outfile, "w") as fasta_out:
        for protein_id, sequence in all_seqs.items():
            if protein_id in ids:
                fasta_out.write(f">{protein_id}\n")
                fasta_out.write(f"{sequence.seq}\n")

    print("fasta writing complete!")
    return


def get_test_train_seqs(
    sp,
    filtered_sp=None,
    train_counts=TEST_TRAIN_COUNTS,
    random_state=CONFIG["random_state"],
    sequences=SEQPATH,
):
    # Split Swiss-Prot into train and test
    train_source = filtered_sp if isinstance(filtered_sp, pd.DataFrame) else sp

    # get unique clusters from the training source and sample
    unique_clusters = train_source["cluster"].unique()
    train_cluster_count = min(len(unique_clusters), train_counts["train"])
    train_clusters = pd.Series(unique_clusters).sample(
        n=train_cluster_count, random_state=random_state
    )

    # Get all rows that belong to the sampled clusters from training source
    cluster_data = train_source[train_source["cluster"].isin(train_clusters)]

    # Sample training instances from these clusters
    sample_size = min(len(cluster_data), train_counts["train"])
    sp_train = cluster_data.sample(n=sample_size, random_state=random_state)

    # For test data, use test_source if provided, otherwise use full sp
    # Don't use any training instances in test data
    if filtered_sp is not None and isinstance(filtered_sp, pd.DataFrame):
        test_sp = sp[~sp.index.isin(sp_train.index)]
    else:
        test_sp = sp[~sp.index.isin(sp_train.index)]

    # split remaining Swiss-Prot into test sets
    sp_test98 = test_sp.sample(n=train_counts["sp_test98"], random_state=random_state)
    sp_remaining = test_sp.drop(sp_test98.index)
    sp_test50 = sp_remaining.sample(n=train_counts["test50"], random_state=random_state)

    training_data = sp_train

    testing_data = {"Test_2%": sp_test98, "Test_50%": sp_test50}

    train_ids = sp_train["uniprot_id"].tolist()
    select_seqs(train_ids, sequences, f"{SP_PATH}/train_swissprot.fasta")

    for test_name, test_df in testing_data.items():
        test_ids = test_df["uniprot_id"].tolist()
        fasta_name = f"{SP_PATH}/{test_name.lower().strip('%')}_swissprot.fasta"
        select_seqs(test_ids, sequences, fasta_name)

    return training_data, testing_data


def plot_length_histogram(df, prefix):
    plt.figure(figsize=(5, 5))
    plt.hist(df["length"])
    plt.xlabel("Length")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"{prefix}length_histogram.png")


swissprots = pd.read_csv(CONFIG["swiss_prot_res"], sep="\t")

swissprots = swissprots[swissprots["is_bacteria"] == "Yes"]
swissprots = swissprots[swissprots["AntiFam count"] == 0]
swissprots["source"] = "Swiss-Prot"

below_100 = swissprots[swissprots["length"] <= 100]

training_data, testing_data = get_test_train_seqs(swissprots, below_100)
os.makedirs(
    f"{SP_PATH}/plots",
    exist_ok=True,
)
plot_length_histogram(training_data, f"{SP_PATH}/plots/Train_")
plot_length_histogram(testing_data["Test_2%"], f"{SP_PATH}/plots//Test2_")
plot_length_histogram(testing_data["Test_50%"], f"{SP_PATH}/plots//Test50_")


def pfam_coverage(df):
    print("Pfam coverage: ", len(df[df["Pfam count"] > 0]) / len(df) * 100)


pfam_coverage(training_data)
pfam_coverage(testing_data["Test_2%"])
pfam_coverage(testing_data["Test_50%"])


def save_df(df, prefix):
    df.to_csv(f"{prefix}info.tsv", sep="\t")


os.makedirs(
    f"{SP_PATH}/data",
    exist_ok=True,
)
save_df(training_data, f"{SP_PATH}/data/Train_")
save_df(testing_data["Test_2%"], f"{SP_PATH}/data/Test2_")
save_df(testing_data["Test_50%"], f"{SP_PATH}/data/Test50_")
