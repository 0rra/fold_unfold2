import argparse
import json
import os
import re
import statistics
from pathlib import Path

import numpy as np
import pandas as pd


# parse esmfold/log files
# get sequence length, pLDDT, pTM
def open_log(dir_name):
    with open(f"{dir_name}/log.txt") as fh:
        log_match = re.compile(
            r"for\s+(\S+)\s+with\s+length\s+(\d+),\s+pLDDT\s+([\d.]+),\s+pTM\s+([\d.]+)"
        )
        scores = []
        for line in fh:
            if log_match.search(line):
                info = line.split("|")[-1]
                res = log_match.search(info)
                res.groups()

                length = res.groups()[1]
                pLDDT = res.groups()[2]
                pTM = res.groups()[3]
                scores = [length, pLDDT, pTM]
            else:
                continue
    if not scores:
        print(f"warning: regex match fail for {dir_name}")
    return scores


def read_logs(path):
    file_path = Path(path)
    seqs = {}
    for f in sorted(file_path.iterdir()):
        if f.is_dir():
            # Extract sequence name by removing 'esmfold_' prefix
            seq_name = f.name.replace("esmfold_", "")
            scores = open_log(f)
            if scores:
                seqs[seq_name] = scores
    return seqs


def write_esm_res(log_dir, res_csv):
    seq_info = read_logs(log_dir)
    df_esm = pd.DataFrame.from_dict(
        seq_info, orient="index", columns=["length", "pLDDT mean", "pTM"]
    )
    df_esm.to_csv(res_csv)


### read colabfold results


def compute_d0(n_res):
    clipped_n_res = max(n_res, 19)  # Clip to greater equals 19
    d0 = 1.24 * (clipped_n_res - 15) ** (1 / 3) - 1.8
    return d0


def calc_ptm(pae_matrix, n_res):
    d0 = compute_d0(n_res)
    conf_matrix = 1 / (1 + (pae_matrix / d0) ** 2)
    row_averages = np.mean(conf_matrix, axis=1)
    ptm = np.max(row_averages)
    return ptm


def read_scores(base_path):
    # read scores from all colabfold subdirectories
    file_path = Path(base_path)
    seqs = {}
    res_plddt = {}
    res_pae = {}

    # Find all subdirectories in the colabfold directory
    subdirs = [d for d in file_path.iterdir() if d.is_dir()]

    for subdir in subdirs:
        # Look for score files in each subdirectory
        score_files = list(subdir.glob("*scores_rank_001_alphafold2_ptm_model*.json"))

        for f in score_files:
            # Extract sequence name
            seq = f.name.split("_scores")[0]
            scores, plddt_per_res, pae_per_res = read_json(str(f))
            seqs[seq] = scores
            res_plddt[seq] = plddt_per_res
            res_pae[seq] = pae_per_res

    return seqs, res_plddt, res_pae


def read_json(json_path):
    print(f"Reading: {json_path}")
    with open(json_path) as json_data:
        pae_j = json.load(json_data)
        average_plddt = statistics.mean(pae_j["plddt"])
        length = len(pae_j["plddt"])
        ptm = pae_j["ptm"]
        plddt_per_res = pae_j["plddt"]
        pae_per_res = pae_j["pae"]
        max_pae = pae_j["max_pae"]
        average_pae = np.mean(pae_j["pae"])
        alt_ptm = calc_ptm(np.array(pae_j["pae"]), length)
        scores = length, average_plddt, ptm, max_pae, average_pae, alt_ptm
    return scores, plddt_per_res, pae_per_res


def write_colabfold_res(res_dir, res_csv, plddt_json, pae_json):
    seqs, res_plddt, res_pae = read_scores(res_dir)

    if not seqs:
        print(f"Warning: No colabfold results found in {res_dir}")
        return

    df_colabfold = pd.DataFrame.from_dict(
        seqs,
        orient="index",
        columns=["length", "pLDDT mean", "pTM", "max pae", "mean pae", "pae pTM"],
    )
    df_colabfold.to_csv(res_csv)
    with open(plddt_json, "w") as f:
        f.write(json.dumps(res_plddt))
    with open(pae_json, "w") as f:
        f.write(json.dumps(res_pae))


### read af3 results
def read_res(res_dir):
    seq_info = {}
    sub_dirs = [f.path for f in os.scandir(res_dir) if f.is_dir()]
    print(f"Processing AF3 subdirectories: {sub_dirs}")

    for seq_dir in sub_dirs:
        seq_name = os.path.basename(seq_dir).upper()

        # Build paths with proper case handling
        conf_path = f"{seq_dir}/{seq_name.lower()}_confidences.json"
        summary_path = f"{seq_dir}/{seq_name.lower()}_summary_confidences.json"
        data_path = f"{seq_dir}/{seq_name.lower()}_data.json"

        # Check if files exist
        if not os.path.exists(conf_path):
            print(f"Warning: Missing confidences file for {seq_name}")
            continue

        with open(conf_path) as cj:
            confidences = json.load(cj)
            atom_plddt = confidences["atom_plddts"]
            fixed_seq_name = seq_name.replace("SEQ", "seq")
            seq_info[fixed_seq_name] = {"pLDDT": atom_plddt}
            seq_info[fixed_seq_name]["pLDDT mean"] = sum(atom_plddt) / len(atom_plddt)
            seq_info[fixed_seq_name]["pae"] = confidences["pae"]
            seq_info[fixed_seq_name]["max pae"] = np.max(confidences["pae"])
            seq_info[fixed_seq_name]["mean pae"] = np.mean(confidences["pae"])
            seq_info[fixed_seq_name]["pae pTM"] = calc_ptm(
                np.array(confidences["pae"]), len(atom_plddt)
            )

        if os.path.exists(summary_path):
            with open(summary_path) as sj:
                summary = json.load(sj)
                seq_info[fixed_seq_name]["pTM"] = summary["ptm"]
        else:
            print(f"Warning: Missing summary file for {seq_name}")
            seq_info[fixed_seq_name]["pTM"] = None

        if os.path.exists(data_path):
            with open(data_path) as dj:
                data = json.load(dj)
                for seq in data["sequences"]:
                    length = len(seq["protein"]["sequence"])
                seq_info[fixed_seq_name]["length"] = length
        else:
            print(f"Warning: Missing data file for {seq_name}")
            seq_info[fixed_seq_name]["length"] = len(atom_plddt)

    return seq_info


def write_af3_res(res_dir, csv_file, plddt_json, pae_json):
    # process all af3 batch directories
    batch_dirs = [f.path for f in os.scandir(res_dir) if f.is_dir()]
    print(f"Found AF3 batch directories: {batch_dirs}")

    af3_res = {}
    for d in batch_dirs:
        new_seq = read_res(d)
        af3_res.update(new_seq)

    if not af3_res:
        print(f"Warning: No AF3 results found in {res_dir}")
        return

    df_af3 = pd.DataFrame.from_dict(
        af3_res,
        orient="index",
        columns=["length", "pLDDT mean", "pTM", "max pae", "mean pae", "pae pTM"],
    )
    df_af3.to_csv(csv_file)

    atom_plddts = {}
    atom_paes = {}
    for seq in af3_res:
        atom_plddts[seq] = af3_res[seq]["pLDDT"]
        atom_paes[seq] = af3_res[seq]["pae"]
    with open(plddt_json, "w") as f:
        f.write(json.dumps(atom_plddts))
    with open(pae_json, "w") as f:
        f.write(json.dumps(atom_paes))


def parse_structures(input_res_dir, parse_res_dir, prefix):
    Path(parse_res_dir).mkdir(parents=True, exist_ok=True)

    # Parse ESMFold results
    esmfold_dir = Path(input_res_dir) / "esmfold"
    if esmfold_dir.exists():
        print(f"\nProcessing ESMFold results from {esmfold_dir}")
        write_esm_res(str(esmfold_dir), f"{parse_res_dir}/{prefix}esmfold_info.csv")
    else:
        print(f"Warning: ESMFold directory not found at {esmfold_dir}")

    # Parse AF3 results
    af3_dir = Path(input_res_dir) / "af3"
    if af3_dir.exists():
        print(f"\nProcessing AF3 results from {af3_dir}")
        write_af3_res(
            str(af3_dir),
            f"{parse_res_dir}/{prefix}af3_info.csv",
            f"{parse_res_dir}/{prefix}af3_plddts.json",
            f"{parse_res_dir}/{prefix}af3_pae.json",
        )
    else:
        print(f"Warning: AF3 directory not found at {af3_dir}")

    # Parse ColabFold results
    colabfold_dir = Path(input_res_dir) / "colabfold"
    if colabfold_dir.exists():
        print(f"\nProcessing ColabFold results from {colabfold_dir}")
        write_colabfold_res(
            str(colabfold_dir),
            f"{parse_res_dir}/{prefix}colabfold_info.csv",
            f"{parse_res_dir}/{prefix}colabfold_plddts.json",
            f"{parse_res_dir}/{prefix}colabfold_pae.json",
        )
    else:
        print(f"Warning: ColabFold directory not found at {colabfold_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse structure prediction results from ESMFold, AF3, and ColabFold"
    )
    parser.add_argument(
        "structure_dir",
        help="Path to directory containing structure predictions (with af3/colabfold/esmfold subdirs)",
    )
    parser.add_argument("parsed_dir", help="Path to directory for parsed results")
    parser.add_argument(
        "--prefix",
        dest="results_prefix",
        help="Desired prefix for parsed result files",
        default="",
    )
    args = parser.parse_args()

    print(f"Structure directory: {args.structure_dir}")
    print(f"Output directory: {args.parsed_dir}")
    print(f"Results prefix: {args.results_prefix}")

    parse_structures(
        args.structure_dir,
        args.parsed_dir,
        args.results_prefix,
    )

    print("\nParsing complete!")


if __name__ == "__main__":
    main()
