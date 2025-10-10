import os
import csv
from pathlib import Path


def extract_folder_id(folder_name):
    """Extract folder_id like NN083_2 from HE_NN083_2_"""
    parts = folder_name.strip("_").split("_")
    return ".".join(parts[1:3])

def sort_key(row):
    """Extract numeric value from folder_id if possible, else use string"""
    folder_id = row[0]
    try:
        return int(''.join([c for c in folder_id if c.isdigit()]) or 0)
    except ValueError:
        return folder_id

def new_csv(base_path: str | Path, output_path: str | Path, cell_type: str):
    
    all_rows = []
    header = ["Biopsi", f"nbr_pixels_{cell_type}_cells", "nbr_tissue_pixels_in_mask", f"nbr_nuclei_{cell_type}_cells", "tot_nuclei_wsi", f"nbr_{cell_type}_cells_per_mm2", "tot_cells_per_mm2", f"percent_pixels_{cell_type}_cells", f"percent_nuclei_{cell_type}_cells"]
    new_header = ["Biopsi", f"percent_pixels_{cell_type}_cells", f"percent_nuclei_{cell_type}_cells", f"nbr_nuclei_{cell_type}_cells", f"nbr_{cell_type}_cells_per_mm2", "tot_nuclei_wsi", "tot_cells_per_mm2", f"nbr_pixels_{cell_type}_cells", "nbr_tissue_pixels_in_mask"]

    # Loop through all folders in base_path
    for folder_name in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder_name)
        
        if os.path.isdir(folder_path):
            stats_path = os.path.join(folder_path, f"stats_{cell_type}", "stats.csv")
            if not os.path.exists(stats_path):
                continue
            
            with open(stats_path, "r", newline="") as f:
                reader = list(csv.reader(f))
                
                # Collect second row only
                if len(reader) > 1:
                    folder_id = extract_folder_id(folder_name)
                    row = [folder_id] + reader[1]
                    all_rows.append(row)

    # Sort rows by folder_id (smallest number first)
    all_rows.sort(key=sort_key)

    # Build a mapping from old index to new order
    index_map = [header.index(col) for col in new_header]

    # Reorder each row
    reordered_rows = []
    for row in all_rows:
        # Pad row if it's shorter than expected
        row = row + [""] * (len(header) - len(row))
        reordered_rows.append([row[i] for i in index_map])

    # Write to csv file
    if new_header:
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(new_header)
            writer.writerows(reordered_rows)

    print(f"Collected {len(reordered_rows)} rows into {output_path}, sorted by folder_id")

if __name__ == "__main__":
    # cells= ["neutrophil", "lymphocyte", "plasma", "eosinophil", "connective", "epithelial"]
    cell_type ="eosinophil"
    base_path = "/media/jenny/Expansion/MM_HE_results"
    output_path = f"/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_{cell_type}_HE.csv"
    new_csv(base_path, output_path, cell_type)
