import argparse

import pandas as pd


def annotate_cluster(cluster_file, df):
    cluster_map = {}
    with open(cluster_file) as cf:
        for line in cf:
            cluster, seqid = line.strip().split("\t")
            cluster_map[seqid] = cluster
    df["cluster"] = df["uniprot_id"].map(cluster_map)
    # map cluster
    print(df.head())
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Parse mmseqs results file and add to existing CSV/JSON"
    )
    parser.add_argument("input_mmseqs", help="Input mmseqs results file")
    parser.add_argument("input_data", help="Input CSV or JSON file to update")
    parser.add_argument("output_file", help="Output file (will match input format)")
    args = parser.parse_args()

    FileIn = args.input_mmseqs
    InputData = args.input_data
    OutputFile = args.output_file

    df = pd.read_csv(InputData, sep="\t")
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "id"})

    df = annotate_cluster(FileIn, df)
    df.to_csv(OutputFile, sep="\t", index=False)


if __name__ == "__main__":
    main()
