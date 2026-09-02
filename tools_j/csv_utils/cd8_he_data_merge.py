import csv
import os
import re

# Paths to CSV files
file1_path = "/media/jenny/Expansion/CellPopEstimatesIHC_QuPathAnalyses.csv"
file2_path = "/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_HE.csv"

# Output
output_path = "/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_cd8_he.csv"

decimals_per_column = {
    "Biopsi": None,     
    "CD8PosPixelsPercent_WholeTissue": int(4),
    "CD8PosCellsPercent_WholeTissue": int(3),
    "NumberOfCD8PosCells_WholeTissue": int(0), 
    "NumberOfCD8PosCellsPermm2_WholeTissue" : int(2),
    "NumberOfCells_CD8section_WholeTissue" : int(0),
    "NumberOfCellsPermm2_CD8section_WholeTissue" : int(3),
}

file2_to_file1 = {
    "Biopsi": "Biopsi",           
    "CD8PosPixelsPercent_WholeTissue" : "percent_pixels_lymphocyte_cells", 
    "CD8PosCellsPercent_WholeTissue": "percent_nuclei_lymphocyte_cells", 
    "NumberOfCD8PosCells_WholeTissue" : "nbr_nuclei_lymphocyte_cells", 
    "NumberOfCD8PosCellsPermm2_WholeTissue" : "nbr_lymphocyte_cells_per_mm2", 
    "NumberOfCells_CD8section_WholeTissue" : "tot_nuclei_wsi", 
    "NumberOfCellsPermm2_CD8section_WholeTissue" : "tot_cells_per_mm2"
}

def extract_columns(row, col_names):
    """Extract only the selected columns from a row"""
    return {col: row.get(col, "") for col in col_names}


def format_value(value, col_name):
    """Format numbers based on column-specific decimals"""
    try:
        num = float(value)
        decimals = decimals_per_column.get(col_name, None)
        if decimals is None:
            return value  # leave unchanged
        return f"{num:.{decimals}f}"
    except (ValueError, TypeError):
        return value  # leave unchanged


# Read headers from file1
with open(file1_path, newline="") as f1:
    reader1 = csv.reader(f1)
    headers1 = next(reader1)
    Biopsi_index = headers1.index("Biopsi")
    selected_columns = ["Biopsi"] + headers1[Biopsi_index+1 : Biopsi_index+7]

# Prepare output columns: source + selected_columns
output_columns = ["source"] + selected_columns

data = []

# Function to process a file
def process_file(file_path, source_label, column_map=None):
    with open(file_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            aligned_data = {}
            for file1_col in selected_columns:
                # Map file2 column if needed
                col = column_map.get(file1_col, file1_col) if column_map else file1_col
                value = row.get(col, "")
                aligned_data[file1_col] = format_value(value, file1_col)
            aligned_data["source"] = source_label
            # Helpers for sorting
            aligned_data["_Biopsi_sort"] = aligned_data["Biopsi"]
            aligned_data["_file_index"] = 0 if source_label == "CD8" else 1
            data.append(aligned_data)

# Process both files
process_file(file1_path, "CD8")
process_file(file2_path, "HE", file2_to_file1)
print("Before sort:", len(data))
# Sort by Biopsi, then by file index
data.sort(key=lambda x: (x["_Biopsi_sort"], x["_file_index"]))
print("After sort:", len(data))
# Write output CSV
os.makedirs(os.path.dirname(output_path), exist_ok=True)
with open(output_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=output_columns)
    writer.writeheader()
    for row in data:
        writer.writerow({k: v for k, v in row.items() if k in output_columns})

print(f"Combined {len(data)} rows into {output_path}")
