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


def parse_pdb_plddt(pdb_path):
    """Returns dict of {res_num: plddt} for CA atoms only."""
    plddts = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 11:
                continue
            if parts[2] != "CA":
                continue
            try:
                res_num = int(parts[5])
                b = float(parts[10])
                plddts[res_num] = b
            except (ValueError, IndexError):
                continue
    return plddts


def read_esmfold_results(path):
    base = Path(path)
    seqs = {}
    per_res_plddt = {}

    for d in sorted(base.iterdir()):
        if not d.is_dir():
            continue

        seq_name = d.name.replace("esmfold_", "")

        log_file = d / "log.txt"
        if not log_file.exists():
            print(f"Warning: missing log.txt in {d}")
            continue

        length, plddt_mean, ptm = open_log(str(d))
        length = int(length)
        plddt_mean = float(plddt_mean)
        ptm = float(ptm)

        pdb_files = list(d.glob("*.pdb"))
        if not pdb_files:
            print(f"Warning: no PDB file found in {d}")
            continue

        plddt_dict = parse_pdb_plddt(pdb_files[0])

        # Build full-length array, NaN where residues are missing
        plddt_per_res = [plddt_dict.get(i, float("nan")) for i in range(1, length + 1)]

        # Report any gaps
        missing = [i for i in range(1, length + 1) if i not in plddt_dict]
        if missing:
            print(
                f"Note: {seq_name} missing pLDDT for residue(s) {missing} (likely X or unmodelled)"
            )

        seqs[seq_name] = {
            "length": length,
            "pLDDT mean": plddt_mean,
            "pTM": ptm,
            "max pae": None,
            "mean pae": None,
            "pae pTM": None,
        }
        per_res_plddt[seq_name] = plddt_per_res

    return seqs, per_res_plddt


def write_esm_res(log_dir, res_csv, plddt_json):
    seqs, per_res_plddt = read_esmfold_results(log_dir)

    if not seqs:
        print(f"Warning: No ESMFold results found in {log_dir}")
        return

    df_esm = pd.DataFrame.from_dict(
        seqs,
        orient="index",
        columns=["length", "pLDDT mean", "pTM", "max pae", "mean pae", "pae pTM"],
    )
    df_esm.to_csv(res_csv)

    with open(plddt_json, "w") as f:
        f.write(json.dumps(per_res_plddt))


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


def parse_cif_plddt(cif_path):
    """Extract per-residue pLDDT from CA atoms in an AF3 mmCIF file."""
    plddt_dict = {}
    with open(cif_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            try:
                if parts[3] != "CA":
                    continue
                res_num = int(parts[8])  # label_seq_id (1-indexed residue number)
                b = float(parts[14])  # B-factor = pLDDT
                plddt_dict[res_num] = b
            except (ValueError, IndexError):
                continue
    return plddt_dict


def read_res(res_dir):
    seq_info = {}
    sub_dirs = [f.path for f in os.scandir(res_dir) if f.is_dir()]
    print(f"Processing AF3 subdirectories: {sub_dirs}")

    for seq_dir in sub_dirs:
        seq_name = os.path.basename(seq_dir).upper()

        conf_path = f"{seq_dir}/{seq_name.lower()}_confidences.json"
        summary_path = f"{seq_dir}/{seq_name.lower()}_summary_confidences.json"
        data_path = f"{seq_dir}/{seq_name.lower()}_data.json"

        if not os.path.exists(conf_path):
            print(f"Warning: Missing confidences file for {seq_name}")
            continue

        with open(conf_path) as cj:
            confidences = json.load(cj)

        atom_plddt = confidences["atom_plddts"]
        pae = confidences["pae"]

        # --- Load true sequence length ---
        if os.path.exists(data_path):
            with open(data_path) as dj:
                data = json.load(dj)
                true_len = sum(
                    len(seq["protein"]["sequence"]) for seq in data["sequences"]
                )
        else:
            print(f"Warning: Missing data file for {seq_name}")
            true_len = None

        # --- Convert per-atom pLDDT → per-residue pLDDT ---
        if "atom_residue_index" in confidences:
            # Old AF3 format
            residue_plddt = {}
            for plddt, res_idx, atom_name in zip(
                atom_plddt,
                confidences["atom_residue_index"],
                confidences["atom_atom_name"],
            ):
                if atom_name == "CA":
                    residue_plddt[res_idx] = plddt
            per_res_plddt = [residue_plddt[i] for i in sorted(residue_plddt.keys())]

        elif "token_res_ids" in confidences:
            # New AF3 format — use CA atoms from CIF file
            cif_files = list(Path(seq_dir).glob("*.cif"))
            if not cif_files:
                print(f"Warning: No CIF file found for {seq_name}")
                continue
            n_tokens = len(confidences["token_res_ids"])
            plddt_dict = parse_cif_plddt(cif_files[0])
            per_res_plddt = [
                plddt_dict.get(i, float("nan")) for i in range(1, n_tokens + 1)
            ]
            missing = [i for i in range(1, n_tokens + 1) if i not in plddt_dict]
            if missing:
                print(f"Note: {seq_name} missing pLDDT for residues {missing}")

        else:
            print(f"Warning: Unknown confidence format for {seq_name}")
            continue

        # --- Check for mismatches ---
        if true_len is not None and len(per_res_plddt) != true_len:
            print(
                f"Warning: AF3 length mismatch for {seq_name}: "
                f"data.json={true_len}, CA_count={len(per_res_plddt)}, "
                f"atom_count={len(atom_plddt)}"
            )

        fixed_seq_name = seq_name.replace("SEQ", "seq")
        seq_info[fixed_seq_name] = {
            "pLDDT": per_res_plddt,
            "pLDDT mean": float(np.nanmean(per_res_plddt)),
            "pae": pae,
            "max pae": float(np.max(pae)),
            "mean pae": float(np.mean(pae)),
            "pae pTM": calc_ptm(np.array(pae), len(per_res_plddt)),
            "length": len(per_res_plddt),
        }

        if os.path.exists(summary_path):
            with open(summary_path) as sj:
                summary = json.load(sj)
                seq_info[fixed_seq_name]["pTM"] = summary["ptm"]
        else:
            seq_info[fixed_seq_name]["pTM"] = None

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
        write_esm_res(
            str(esmfold_dir),
            f"{parse_res_dir}/{prefix}esmfold_info.csv",
            f"{parse_res_dir}/{prefix}esmfold_plddts.json",
        )
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
