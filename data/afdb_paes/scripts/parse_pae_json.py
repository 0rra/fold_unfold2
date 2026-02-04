import json
import argparse
import glob
import os
import re
import numpy as np

def compute_d0(n_res):
    clipped_n_res = max(n_res, 19) # Clip to ≥19
    d0 = 1.24 * (clipped_n_res - 15) ** (1/3) - 1.8
    return d0

def calc_ptm(pae_matrix, n_res):
    d0 = compute_d0(n_res)
    conf_matrix = 1 / (1 + (pae_matrix / d0) ** 2)
    row_averages = np.mean(conf_matrix, axis=1)
    ptm = np.max(row_averages)
    return ptm
    

def read_json(json_path):
    filename = os.path.basename(json_path)
    match = re.search(r'AF-(.+?)-predicted_aligned_error_v6\.json', filename)
    if match:
        afdb_id = "AF-" + match.group(1)
    
    with open(json_path) as json_data:
        data = json.load(json_data)
        
        if isinstance(data, list) and len(data) > 0:
            pae_data = data[0]
        else:
            pae_data = data

        pae_matrix = np.array(pae_data["predicted_aligned_error"])
        max_pae = pae_data["max_predicted_aligned_error"]
        
        average_pae = np.mean(pae_matrix)
        
        length = len(pae_matrix)
        
        alt_ptm = calc_ptm(pae_matrix, length)
    
    return afdb_id, max_pae, average_pae, alt_ptm


def process_directory(directory_path, output_file="alphafold_pae_results.tsv"):
    pattern = os.path.join(directory_path, "AF-*-predicted_aligned_error_v6.json")
    json_files = glob.glob(pattern)
    
    if not json_files:
        print(f"No files matching pattern found in {directory_path}")
        return

    
    print(f"Found {len(json_files)} files to process")
    
    results = []

    for json_file in json_files:
        try:
            afdb_id, max_pae, average_pae, alt_ptm = read_json(json_file)
            results.append([afdb_id, max_pae, average_pae, alt_ptm])
#            print(f"Processed: {afdb_id}")
        except Exception as e:
            print(f"Error processing {json_file}: {e}")
    
    with open(output_file, 'w') as tsv_file:
        tsv_file.write("afdb_id\tmax_pae\taverage_pae\talt_ptm\n")
        
        for row in results:
            tsv_file.write(f"{row[0]}\t{row[1]:.4f}\t{row[2]:.4f}\t{row[3]:.4f}\n")
    
    print(f"Results written to: {output_file}")
    print(f"Processed {len(results)} files successfully")


def main():    
    parser = argparse.ArgumentParser(description="Process AlphaFold PAE JSON files")
    parser.add_argument("directory", help="Directory containing AF-*-predicted_aligned_error_v4.json files")
    parser.add_argument("-o", "--output", default="alphafold_pae_results.tsv", 
                       help="Output TSV filename (default: alphafold_pae_results.tsv)")
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a valid directory")
        return
  
    process_directory(args.directory, args.output)


if __name__ == "__main__":
    main()
