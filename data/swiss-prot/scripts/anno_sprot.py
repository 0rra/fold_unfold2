import argparse
import time

import pandas as pd
import requests
from Bio import SeqIO

TAX_FILE = "data/taxonomy/taxonomy_ancestor_2_2025_10_17.tsv"


def get_bact_ids(tax_file):
    bac_taxids = set()
    with open(tax_file) as f:
        for line in f:
            if line.startswith("Taxon"):
                continue
            else:
                line_parts = line.split("\t")
                taxid = line_parts[0]
                lineage = line_parts[1]
                lineage = lineage.split(", ")
                if len(lineage) > 1 and "Bacteria" in lineage[1]:
                    bac_taxids.add(taxid)

    return bac_taxids


# parse out uniprot_id, tax_id, sequence length
def make_seq_df(fasta_file):
    seq_data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        sequence = record.seq
        length = len(sequence)
        uniprot_id = ""
        name = ""

        # extract Uniprot ID and name
        if "|" in header:
            header_parts = header.split("|")
            if len(header_parts) >= 3:
                if header_parts[0] == "sp":
                    uniprot_id = header_parts[1]
                name_parts = header_parts[2].split(" ", 1)
                if len(name_parts) > 1:
                    name_end = name_parts[1]
                    name = name_end.split(" OS=")[0]

        # extract organism tax ID
        ox = "Unknown"
        if "OX=" in header:
            ox = header.split("OX=")[1].split()[0]

        seq_data.append(
            {
                "uniprot_id": uniprot_id,
                "taxid": ox,
                "length": length,
                "name": name,
                # "sequence": str(sequence),
            }
        )

    seq_df = pd.DataFrame(seq_data)
    return seq_df


def anno_tax(seq_df, taxid_list):
    seq_df["is_bacteria"] = seq_df["taxid"].apply(
        lambda x: "Yes" if str(x) in taxid_list else "No"
    )
    return seq_df


def filter_fasta(seq_ids, in_fasta, out_fasta):
    filtered_records = []

    seq_id_set = set(seq_ids)

    for record in SeqIO.parse(in_fasta, "fasta"):
        header = record.description

        uniprot_id = ""
        if "|" in header:
            header_parts = header.split("|")
            if len(header_parts) >= 2:
                uniprot_id = header_parts[1]

        if uniprot_id in seq_id_set:
            filtered_records.append(record)

    SeqIO.write(filtered_records, out_fasta, "fasta")

    return len(filtered_records)


def get_reasons(id):
    url = f"https://rest.uniprot.org/uniprotkb/{id}.json"
    response = requests.get(url)
    annotation_score = None
    protein_existence = None
    caution = None
    is_fragment = "No"

    if response.status_code != 200:
        print(f"Error occured retrieving data for {id}: {response.status_code}")
        return annotation_score, protein_existence, caution, is_fragment

    try:
        data = (
            response.json()
        )
        annotation_score = data.get("annotationScore")
        protein_existence = data.get("proteinExistence")

        if "comments" in data:
            comments = data["comments"]
            for comment in comments:
                if "CAUTION" in comment.get("commentType", "").upper():
                    if "sequenceCautionType" in comment:
                        caution = comment["sequenceCautionType"]
                    elif (
                        "texts" in comment and len(comment["texts"]) > 0
                    ):
                        caution = comment["texts"][0]["value"]

        protein_description = data.get("proteinDescription", {})
        if "fragment" in protein_description.get("flag", "").lower():
            is_fragment = "Yes"

    except Exception as e:
        print(f"Error occured retrieving data for {id}: {e}")
    time.sleep(0.5)
    return annotation_score, protein_existence, caution, is_fragment


def anno_reasons(results_df):
    results_df[
        ["annotation_score", "protein_existence_score", "caution", "is_fragment"]
    ] = results_df["uniprot_id"].apply(lambda acc: pd.Series(get_reasons(acc)))
    return results_df


def main():
    parser = argparse.ArgumentParser(description="Annotate bacteria sequences.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("out_prefix", help="output name for fasta file and output tsv")
    parser.add_argument("--tax_file", default=TAX_FILE)
    args = parser.parse_args()

    input_fasta = args.input_fasta
    tax_file = args.tax_file
    out_prefix = args.out_prefix
    print("getting bacteria ids")
    bact_ids = get_bact_ids(tax_file)

    seq_df = make_seq_df(input_fasta)
    print("parsed sequence fasta to make dataframe")
    seq_df = anno_tax(seq_df, bact_ids)
    print("annotated taxonomy")
    # Save seq_df as tsv
    print("saving without uniprot annotations")
    seq_df.to_csv(f"{out_prefix}_tmp.tsv", sep="\t", index=False)

    # Filter for bacterial sequences
    bact_only = seq_df[seq_df["is_bacteria"] == "Yes"]["uniprot_id"].tolist()
    print("making bact only fasta")
    if bact_only:
        num_filtered = filter_fasta(bact_only, input_fasta, f"{out_prefix}_bact.fasta")
        print(f"Filtered {num_filtered} bacterial sequences to {out_prefix}_bact.fasta")
    else:
        print("No bacterial sequences found.")
    print("getting uniprot annotations")
    seq_df = anno_reasons(seq_df)
    # Save seq_df as tsv
    seq_df.to_csv(f"{out_prefix}.tsv", sep="\t", index=False)

    bact_df = seq_df[seq_df["is_bacteria"] == "Yes"]
    bact_df.to_csv(f"{out_prefix}_bacteria.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
