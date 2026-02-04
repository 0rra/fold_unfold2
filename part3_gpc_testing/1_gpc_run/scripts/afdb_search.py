import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    precision_recall_curve,
)
from sklearn.preprocessing import StandardScaler


def drop_dup_antifams(df):
    df = df.copy()
    df["anf_prefix"] = df["id"].str.extract(r"(ANF\d+)")
    anf = df[df["anf_prefix"].notna()]
    anf = anf[anf.groupby("anf_prefix")["id"].transform("first") == anf["id"]]
    anf = anf[anf["id"] != "ANF00264_seq_1"]

    return pd.concat([anf, df[df["anf_prefix"].isna()]]).drop(columns="anf_prefix")


seed_af_dir = "part1_refoldunfold/2_structure_prediction/parsed_results"
pseudo_af_dir = "part2_gpc_training/2_structure_prediction/parsed_results"
swissprot_dir = "part2_gpc_training/2_structure_prediction/parsed_results"


CONFIG = {
    "feats": ["pae pTM", "pLDDT mean", "length"],
    "pos_label": "AntiFam",
    "random_state": 42,
    "test_counts": {"af_test2": 100, "af_test50": 300},
}


def load_data(directory, prefix, source):
    file_path = f"{directory}/{prefix}colabfold_info.csv"
    df = pd.read_csv(file_path)

    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "id"})

    df["source"] = source

    if source == "AntiFam":
        df["AntiFam source"] = "Generated" if "pseudo" in prefix else "Seed"
        if not prefix:
            df = drop_dup_antifams(df)

    return df


# load datasets
seed_af = load_data(seed_af_dir, "af_res_", "AntiFam")
pseudo_af = load_data(pseudo_af_dir, "pseudo_antifams_", "AntiFam")
train_sp = load_data(swissprot_dir, "train_", "Swiss-Prot")
test50_sp = load_data(swissprot_dir, "test50_", "Swiss-Prot")
test2_sp = load_data(swissprot_dir, "test2_", "Swiss-Prot")

# prepare training data
training_data = pd.concat([seed_af, train_sp])
train_ids = train_sp["id"].tolist()

# prepare test datasets
af_test2 = pseudo_af.sample(
    n=CONFIG["test_counts"]["af_test2"], random_state=CONFIG["random_state"]
)
af_test50 = pseudo_af.sample(
    n=CONFIG["test_counts"]["af_test50"], random_state=CONFIG["random_state"]
)

test2_combined = pd.concat([af_test2, test2_sp])
test50_combined = pd.concat([af_test50, test50_sp])

test_ids = pd.concat([test2_sp, test50_sp])["id"].tolist()
swissprot_ids = train_ids + test_ids

test_data = {
    "Test_50%": test50_combined[CONFIG["feats"] + ["source"]],
    "Test_2%": test2_combined[CONFIG["feats"] + ["source"]],
}

# store data with ids
test_data_with_ids = {
    "Test_50%": test50_combined,
    "Test_2%": test2_combined,
}

# train model
X_train = training_data[CONFIG["feats"]]
y_train = training_data["source"]

print(f"Training with {len(train_ids)} Swiss-Prot sequences")
print(f"Total Swiss-Prot IDs tracked: {len(swissprot_ids)}")


def train_gpc_model(X_train, y_train, random_state=42):
    my_scaler = StandardScaler()
    X_train_scaled = my_scaler.fit_transform(X_train)

    kernel = RBF(length_scale=1.0)

    gpc = GaussianProcessClassifier(
        kernel=kernel, random_state=random_state, max_iter_predict=1000
    )
    gpc.fit(X_train_scaled, y_train)
    return gpc, my_scaler


gpc, scaler = train_gpc_model(X_train, y_train)
positive_class_index = list(gpc.classes_).index("AntiFam")

# evaluate on test sets
for test_name, test_df in test_data.items():
    print(f"\n{'=' * 50}")
    print(f"Evaluating on {test_name}")
    print("=" * 50)

    X_test = test_df[CONFIG["feats"]]
    y_test = test_df["source"]
    X_test_scaled = scaler.transform(X_test)

    y_probs = gpc.predict_proba(X_test_scaled)[:, positive_class_index]

    # find optimal threshold at 100% precision
    precision, recall, thresholds = precision_recall_curve(
        y_test, y_probs, pos_label=CONFIG["pos_label"]
    )

    opt_threshold = next((t for p, t in zip(precision, thresholds) if p == 1.0), 0.5)

    # make predictions
    y_pred = np.where(y_probs >= opt_threshold, CONFIG["pos_label"], "Swiss-Prot")

    # summarise / evaluate
    print(f"Optimal Threshold: {opt_threshold:.4f}\n")
    print("Confusion Matrix:")
    print(confusion_matrix(y_test, y_pred, labels=["Swiss-Prot", CONFIG["pos_label"]]))
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred, target_names=gpc.classes_))

print("using threshold:", opt_threshold)

colab_pae_ptm = {"scaler": scaler, "model": gpc, "threshold": opt_threshold}


def run_model(model_dict, batch_df, mod_feats, name):
    model_scaler = model_dict["scaler"]
    model_gpc = model_dict["model"]
    model_threshold = model_dict["threshold"]

    X_test = batch_df[mod_feats]
    X_scaled = model_scaler.transform(X_test)

    positive_class_index = list(model_gpc.classes_).index(CONFIG["pos_label"])
    y_probs = model_gpc.predict_proba(X_scaled)[:, positive_class_index]

    y_pred = np.where(y_probs > model_threshold, CONFIG["pos_label"], "SwissProt")

    result_df = batch_df.copy()
    result_df["prediction_score"] = y_probs
    result_df["prediction_label"] = y_pred
    result_df["model_name"] = name

    return result_df


# got trained model now need to load...

print("-------------running on afdb data---------------")

afdb_info = "data/afdb_bact_only/2025_03_afdb_bact_info.tsv"

print("removing swissprot ids used for training from afdb eg...", swissprot_ids[3:10])
batches = []

total_rows = sum(1 for _ in open(afdb_info)) - 1  # -1 for header
chunk_size = 10000
total_batches = (total_rows + chunk_size - 1) // chunk_size
batch_num = 0

for batch in pd.read_csv(
    afdb_info,
    chunksize=10000,
    iterator=True,
    sep="\t",  # remove index_col=False to avoid extra index
):
    batch_num += 1
    batch_clean = batch[~batch["uniprot_id"].isin(swissprot_ids)]
    batch_pred = run_model(colab_pae_ptm, batch_clean, CONFIG["feats"], "AFDB")
    batches.append(batch_pred)
    print(f"Running predictions... batch {batch_num}/{total_batches}")

# combine all at the end
print("joining batches")
predictions_df = pd.concat(batches, ignore_index=True)

# remove any accidental duplicate index column
if "Unnamed: 0" in predictions_df.columns:
    predictions_df = predictions_df.drop(columns=["Unnamed: 0"])
print("saving to file")
predictions_df.to_csv(
    "part3_gpc_testing/results/afdb_predictions.tsv",
    sep="\t",
    index=False,
)
