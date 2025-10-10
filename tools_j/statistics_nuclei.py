import pandas as pd
import os
import tifffile
import xml.etree.ElementTree as ET

from pathlib import Path
from pixel_count import count_black_pixels
from read_csv_files import read_metrics_data

def statistics(mask_path: str | Path, csv_pixel_path: str | Path, csv_nuclei_number_path: str | Path, pixel_size_mm2: float, cells= str):
    
    # Number of black pixels in mask represents number of pixels that are tissue in the image
    nbr_tissue_pixels_in_mask = count_black_pixels(mask_path)
    
    # Total number of pixels for each cell category found in wsi
    tot_pixels_df = read_metrics_data(csv_pixel_path, encoding="ISO-8859-1")
    # Total number of nuclei for each cell category found in wsi
    tot_nuclei_df = read_metrics_data(csv_nuclei_number_path, encoding="ISO-8859-1")
    
    # Extract information for all the immune cells combined
    if cells== "all": 
        tot_pixels_immune_cells = tot_pixels_df["neutrophil"].iloc[0] + tot_pixels_df["lymphocyte"].iloc[0] + tot_pixels_df["plasma"].iloc[0] + tot_pixels_df["eosinophil"].iloc[0]
        print("---------------------------------------------------------------------------------------")
        print(f"Total number of pixels that are immune cells in wsi: {tot_pixels_immune_cells}")
        print("---------------------------------------------------------------------------------------")
        tot_nuclei_immune_cells = tot_nuclei_df["neutrophil"].iloc[0] + tot_nuclei_df["lymphocyte"].iloc[0] + tot_nuclei_df["plasma"].iloc[0] + tot_nuclei_df["eosinophil"].iloc[0]
        print(f"Total number of immune cells in wsi: {tot_nuclei_immune_cells}")
        print("---------------------------------------------------------------------------------------")
        tot_nuclei_wsi = tot_nuclei_df["neutrophil"].iloc[0] + tot_nuclei_df["lymphocyte"].iloc[0] + tot_nuclei_df["plasma"].iloc[0] + tot_nuclei_df["eosinophil"].iloc[0] + tot_nuclei_df["connective"].iloc[0] + tot_nuclei_df["epithelial"].iloc[0]
        print(f"Total number of nuclei in wsi: {tot_nuclei_wsi}")
        print("---------------------------------------------------------------------------------------")
        percent_pixels_immune_cells = (tot_pixels_immune_cells/nbr_tissue_pixels_in_mask) * 100
        print(f"% of pixels categorized as immune cells in tissue region of wsi: {percent_pixels_immune_cells:.4f}")
        print("---------------------------------------------------------------------------------------")
        percent_nuclei_immune_cells = (tot_nuclei_immune_cells/tot_nuclei_wsi) * 100
        print(f"% of nuclei detected in wsi categorized as immune cells: {percent_nuclei_immune_cells:.4f}")
        print("---------------------------------------------------------------------------------------")
        tot_immune_cells_per_mm2 = tot_nuclei_immune_cells/(nbr_tissue_pixels_in_mask*pixel_size_mm2)
        print(f"Total number of immune cells per mm^2 in wsi: {tot_immune_cells_per_mm2}")
        print("---------------------------------------------------------------------------------------")
        tot_cells_per_mm2 = tot_nuclei_wsi/(nbr_tissue_pixels_in_mask*pixel_size_mm2)
        print(f"Total number of cells per mm^2 in wsi: {tot_cells_per_mm2}")
        print("---------------------------------------------------------------------------------------")
        print(f"Total number of tissue pixels in wsi: {nbr_tissue_pixels_in_mask}")
        print("---------------------------------------------------------------------------------------")
        pd.set_option('display.max_rows', None)     # Show all rows
        pd.set_option('display.max_columns', None)  # Show all columns
        pd.set_option('display.width', 1000)        # Adjust the width of the display
        info_df = pd.DataFrame({"tot_pixels_immune_cells": tot_pixels_immune_cells,
                                "nbr_tissue_pixels_in_mask": nbr_tissue_pixels_in_mask,
                                "tot_nuclei_immune_cells":tot_nuclei_immune_cells,
                                "tot_nuclei_wsi": tot_nuclei_wsi,
                                "tot_immune_cells_per_mm2":tot_immune_cells_per_mm2,
                                "tot_cells_per_mm2": tot_cells_per_mm2,
                                "percent_pixels_immune_cells":percent_pixels_immune_cells,
                            "percent_nuclei_immune_cells": percent_nuclei_immune_cells}, index=[0])
        info_df.to_csv(save_stats_path, index=False)
    # Extraxt information for each cell type
    else:
        for cell in cells:
            print(f"Calculates statistics for {cell}")
            tot_pixels_specific_cell = tot_pixels_df[cell].iloc[0]
            print("---------------------------------------------------------------------------------------")
            print(f"Total number of pixels that are {cell} cells in wsi: {tot_pixels_specific_cell}")
            print("---------------------------------------------------------------------------------------")
            tot_nuclei_specific_cell = tot_nuclei_df[cell].iloc[0]
            print(f"Total number of {cell} cells in wsi: {tot_nuclei_specific_cell}")
            print("---------------------------------------------------------------------------------------")
            tot_nuclei_wsi = tot_nuclei_df["neutrophil"].iloc[0] + tot_nuclei_df["lymphocyte"].iloc[0] + tot_nuclei_df["plasma"].iloc[0] + tot_nuclei_df["eosinophil"].iloc[0] + tot_nuclei_df["connective"].iloc[0] + tot_nuclei_df["epithelial"].iloc[0]
            print(f"Total number of nuclei in wsi: {tot_nuclei_wsi}")
            print("---------------------------------------------------------------------------------------")
            percent_pixels_specific_cell = (tot_pixels_specific_cell/nbr_tissue_pixels_in_mask) * 100
            print(f"% of pixels categorized as {cell} cells in tissue region of wsi: {percent_pixels_specific_cell:.4f}")
            print("---------------------------------------------------------------------------------------")
            percent_nuclei_specific_cell = (tot_nuclei_specific_cell/tot_nuclei_wsi) * 100
            print(f"% of nuclei detected in wsi categorized as {cell} cells: {percent_nuclei_specific_cell:.4f}")
            print("---------------------------------------------------------------------------------------")
            tot_specific_cell_per_mm2 = tot_nuclei_specific_cell/(nbr_tissue_pixels_in_mask*pixel_size_mm2)
            print(f"Total number of {cell} cells per mm^2 in wsi: {tot_specific_cell_per_mm2}")
            print("---------------------------------------------------------------------------------------")
            tot_cells_per_mm2 = tot_nuclei_wsi/(nbr_tissue_pixels_in_mask*pixel_size_mm2)
            print(f"Total number of cells per mm^2 in wsi: {tot_cells_per_mm2}")
            print("---------------------------------------------------------------------------------------")
            print(f"Total number of tissue pixels in wsi: {nbr_tissue_pixels_in_mask}")
            print("---------------------------------------------------------------------------------------")
            pd.set_option('display.max_rows', None)     # Show all rows
            pd.set_option('display.max_columns', None)  # Show all columns
            pd.set_option('display.width', 1000)        # Adjust the width of the display
            info_df = pd.DataFrame({f"tot_pixels_{cell}_cells": tot_pixels_specific_cell,
                                    "nbr_tissue_pixels_in_mask": nbr_tissue_pixels_in_mask,
                                    f"tot_nuclei_{cell}_cells":tot_nuclei_specific_cell,
                                    "tot_nuclei_wsi": tot_nuclei_wsi,
                                    f"tot_{cell}_cells_per_mm2":tot_specific_cell_per_mm2,
                                    "tot_cells_per_mm2": tot_cells_per_mm2,
                                    f"percent_pixels_{cell}_cells":percent_pixels_specific_cell,
                                f"percent_nuclei_{cell}_cells": percent_nuclei_specific_cell}, index=[0])
            # Output directory
            save_stats_path = Path(f"/media/jenny/Expansion/MM_HE_results/xxx/stats_{cell}/")
            if not save_stats_path.exists(): 
                save_stats_path.mkdir(parents=True)
                print(f"Directory {save_stats_path} was created")
            save_stats_path = os.path.join(save_stats_path, f"stats.csv")
            info_df.to_csv(save_stats_path, index=False)

def extract_pixel_size(tif_path: str | Path):

    with tifffile.TiffFile(tif_path) as tif:
        omexml = tif.ome_metadata  # Extract OME-XML as string

    # Parse XML
    root = ET.fromstring(omexml)

    # Find Pixels tag 
    pixels = None
    for elem in root.iter():
        if elem.tag.endswith("Pixels"):
            pixels = elem
            break

    if pixels is None:
        raise ValueError("No <Pixels> element found in OME-XML metadata")
    
    # Find pixel size
    physical_size_x = float(pixels.attrib["PhysicalSizeX"])*0.001 # Convert from µm to mm
    physical_size_y = float(pixels.attrib["PhysicalSizeY"])*0.001 # Convert from µm to mm

    pixel_size = physical_size_x*physical_size_y # mm²

    print("PhysicalSizeX:", physical_size_x)
    print("PhysicalSizeY:", physical_size_y)
    return pixel_size


if __name__ == "__main__":
    tif_path = "/media/.../HE_xxx.tif"
    pixel_size = extract_pixel_size(tif_path)
    cells= ["neutrophil", "lymphocyte", "plasma", "eosinophil", "connective", "epithelial"]
    mask_path = "/media/jenny/Expansion/MM_HE_masks/xxx/xxx_full_mask.png"
    csv_nuclei_number_path = "/media/jenny/Expansion/MM_HE_results/xxx/counts/"
    csv_pixel_path = "/media/jenny/Expansion/MM_HE_results/xxx/count_pixels/"
   
    statistics(mask_path, csv_pixel_path, csv_nuclei_number_path, pixel_size, cells)