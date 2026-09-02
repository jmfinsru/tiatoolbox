import pandas as pd
import os
import tifffile
import xml.etree.ElementTree as ET

from pathlib import Path
from pixel_count import count_black_pixels
from read_csv_files import read_metrics_data

def statistics_one_image(mask_path: str | Path, csv_pixel_path: str | Path, csv_nuclei_number_path: str | Path, pixel_size_mm2: float, cells= str):
    
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
            save_stats_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/HE_MIMB755_2B_150425/statistics_cells/stats_{cell}/")
            if not save_stats_path.exists(): 
                save_stats_path.mkdir(parents=True)
                print(f"Directory {save_stats_path} was created")
            save_stats_path = os.path.join(save_stats_path, f"stats.csv")
            info_df.to_csv(save_stats_path, index=False)

def statistics_whole_folder(root_wsi_dir: str | Path, mask_root_dir: str | Path, csv_stats_root_dir: str | Path, cells= str):
    
    # Find .csv files with pixel and nuclei statistics for each biopsy
    csv_stats_root_dir = Path(csv_stats_root_dir)
    files_of_interest = "HE_MIMB77*"
    for collection_dir in sorted(csv_stats_root_dir.glob(files_of_interest)):
        base_name = collection_dir.name
        print(f"Biopsy processed: {base_name}")

        count_pixels_dir = collection_dir / "count_pixels/"
        # csv_pixel_paths = sorted(count_pixels_dir.glob("*.csv"))
        # if len(csv_pixel_paths) != 1:
        #     print(f"The number of csv files are: {len(csv_pixel_paths)}, not 1 as expected.")
        # csv_pixel_path = str(csv_pixel_paths[0])
        
        count_nuclei_dir = collection_dir / "counts/"

        # csv_nuclei_number_paths = sorted(count_nuclei_dir.glob("*.csv"))
        # if len(csv_nuclei_number_paths) != 1:
        #     print(f"The number of csv files are: {len(csv_nuclei_number_paths)}, not 1 as expected.")
        # csv_nuclei_number_path = str(csv_nuclei_number_paths[0])

        mask_root_dir = Path(mask_root_dir)
        mask_collection_dir = sorted(mask_root_dir.glob(f"{base_name}*"))
        if len(mask_collection_dir) != 1:
            print(f"The number of directories found are {len(mask_collection_dir)}, not 1 as expected.")
        mask_collection_dir = mask_collection_dir[0]
        mask_path = sorted(mask_collection_dir.glob("*no_folds_full_mask.tif"))
        mask_path = Path(mask_path[0])
        
        # Find WSI
        wsi_dir = root_wsi_dir / f"{base_name}.vsi.Collection/"
        # Find all tif files under this collection folder whose name starts with the base name
        wsi_path = next(wsi_dir.rglob(f"{base_name}_20x_BF_01.tif"), None)
        
        if wsi_path is None:
            raise FileNotFoundError(f"No matching .tif found in {wsi_dir}")
        pixel_size_mm2 = extract_pixel_size(wsi_path)

        # Number of black pixels in mask represents number of pixels that are tissue in the image
        nbr_tissue_pixels_in_mask = count_black_pixels(mask_path)
        
        # Total number of pixels for each cell category found in wsi
        tot_pixels_df = read_metrics_data(str(count_pixels_dir), encoding="ISO-8859-1")
        # Total number of nuclei for each cell category found in wsi
        tot_nuclei_df = read_metrics_data(str(count_nuclei_dir), encoding="ISO-8859-1")
        
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
                save_stats_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/{base_name}/statistics_cells/stats_{cell}/")
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
    tif_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF/HE_MIMB755_2B_150425.vsi.Collection/HE_MIMB755_2B_150425_Layer3-20x_BF_04/HE_MIMB755_2B_150425_Layer3-20x_BF_04.tif"
    pixel_size = extract_pixel_size(tif_path)
    cells= ["neutrophil", "lymphocyte", "plasma", "eosinophil", "connective", "epithelial"]
    mask_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_masks/HE_MIMB755_2B_150425_Layer3-20x_BF_04/HE_MIMB755_2B_150425_Layer3-20x_BF_04_no_folds_full_mask.tif"
    csv_nuclei_number_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_results/HE_MIMB755_2B_150425/counts/"
    csv_pixel_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_results/HE_MIMB755_2B_150425/count_pixels/"
    
    # root_wsi_dir = Path("/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF/")
    # mask_root_dir = "/media/jenny/Expansion/MetoxyLacc_HE_20x_masks/"
    # csv_stats_root_dir = "/media/jenny/Expansion/MetoxyLacc_HE_20x_results/"
    statistics_one_image(mask_path, csv_pixel_path, csv_nuclei_number_path, pixel_size, cells)
    # statistics_whole_folder(root_wsi_dir=root_wsi_dir, mask_root_dir=mask_root_dir, csv_stats_root_dir=csv_stats_root_dir, cells = cells)