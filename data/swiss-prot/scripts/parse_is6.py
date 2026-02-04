import argparse

import pandas as pd


def anno_is6(is6_path, result_df):
    is6 = pd.read_csv(is6_path, sep="\t", index_col=False, header=None)

    # make cols and split accession col
    header = [
        "UniProt ID",
        "md5",
        "Sequence length",
        "DB",
        "Family ID",
        "Desc",
        "Start",
        "End",
        "Evalue",
        "2",
        "3",
        "InterPro ID",
        "IP name",
        "4",
        "5",
    ]
    # add header
    is6.columns = header

    is6 = is6[["UniProt ID", "DB", "Family ID", "Desc", "Start", "End"]]
    # Check if any values contain "sp|"
    if is6["UniProt ID"].str.contains("sp|", regex=False).any():
        print("fixing fasta header id...")
        is6["UniProt ID"] = is6["UniProt ID"].str.split("|").str[1]
    print(is6)
    ids = result_df[["uniprot_id"]].drop_duplicates()

    # count pfam/antifam per id
    counts = (
        is6.groupby(["UniProt ID", "DB"]).size().unstack(fill_value=0).reset_index()
    )

    # merge back to keep all ids, even if counts are missing
    summary = ids.merge(
        counts, right_on="UniProt ID", left_on="uniprot_id", how="left"
    ).fillna(0)

    # make sure both columns exist
    for col in ["Pfam", "AntiFam", "MobiDB-lite"]:
        if col not in summary.columns:
            summary[col] = 0

    # clean column names
    summary = summary.rename(
        columns={
            "Pfam": "Pfam count",
            "AntiFam": "AntiFam count",
            "MobiDB-lite": "MobiDB count",
        }
    )
    summary = summary.drop(columns=["UniProt ID"])

    result_df = result_df.merge(summary, how="left", on="uniprot_id")

    result_df["IS6 result"] = ~(
        (result_df["Pfam count"] == 0)
        & (result_df["AntiFam count"] == 0)
        & (result_df["MobiDB count"] == 0)
    )
    print(
        f"Checking...number seqs with pfam match: {len(result_df[result_df['Pfam count'] > 0])}\nnumber seqs with antifam match: {len(result_df[result_df['AntiFam count'] > 0])}"
    )
    return is6, result_df


def main():
    parser = argparse.ArgumentParser(
        description="Parse is6 results file and add to existing CSV/JSON"
    )
    parser.add_argument("input_is6", help="Input is6 results file")
    parser.add_argument("input_data", help="Input CSV or JSON file to update")
    parser.add_argument("output_file", help="Output file (will match input format)")
    args = parser.parse_args()

    FileIn = args.input_is6
    InputData = args.input_data
    OutputFile = args.output_file

    df = pd.read_csv(InputData, sep="\t")
    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "id"})

    if "Pfam count" in df.columns:
        df = df.drop(columns=["Pfam count", "AntiFam count", "MobiDB count"])
    print(df.columns)
    print(df.head())
    is6_df, df = anno_is6(FileIn, df)
    df.to_csv(OutputFile, sep="\t", index=False)
    is6_df.to_csv("is6_results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
