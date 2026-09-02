import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd

from scipy import stats
from scipy.stats import wilcoxon, ttest_rel
from pathlib import Path
from itertools import product
from matplotlib.colors import to_rgba

def extract_HE_and_CD8_data(column_list: list, file_path_HE: str | Path, file_path_CD8: str | Path):
    
    # Load data files as a Dataframe
    df_HE = pd.read_csv(file_path_HE)
    df_CD8 = pd.read_csv(file_path_CD8)
    
    # Regex patterns for the Biopsies
    pattern_before = re.compile(r"\.[A-Z]$")
    pattern_after = re.compile(r"\.\d+[A-Z]$")

    results = []
    # Loop through all biopsies in the CD8 Dataframe
    for biopsi in df_CD8["Biopsi"]:
        CD8_dict = {}
        skip_outer = False
        for column in column_list: # Loop through columns of interest to extract the data
            col = df_CD8.loc[df_CD8["Biopsi"]== biopsi, column].dropna().values # Drop NaN values
            if len(col)==0: # Happens when NaN is dropped. skip_outer makes it skip to the next Biopsi
                skip_outer = True 
            else:
                CD8_dict[column] = float(col)
        if skip_outer:
            continue
        
        # Group the biopsies into before and after treatment
        if pattern_before.search(biopsi):
            group = "before"
        elif pattern_after.search(biopsi):
            group = "after"
        else:
            None
        
        # Only include biopsies found in both input files in the results
        if biopsi in df_HE["Biopsi"].values:
            print(biopsi)
            percent = df_HE.loc[df_HE["Biopsi"] == biopsi, f"percent_nuclei_{cell}_cells"].values
            print(percent)
            nbr_type_cells_per_mm2 = df_HE.loc[df_HE["Biopsi"] == biopsi, f"nbr_{cell}_cells_per_mm2"].values
            nbr_nuclei_type_cells = df_HE.loc[df_HE["Biopsi"] == biopsi, f"nbr_nuclei_{cell}_cells"].values
            percent_pixels_type_cells = df_HE.loc[df_HE["Biopsi"] == biopsi, f"percent_pixels_{cell}_cells"].values
            tot_cells_per_mm2 = df_HE.loc[df_HE["Biopsi"] == biopsi, "tot_cells_per_mm2"].values
            tot_nuclei_wsi = df_HE.loc[df_HE["Biopsi"] == biopsi, "tot_nuclei_wsi"].values
            CD8_dict.update({
                    "Biopsi": biopsi,
                    "Group": group,
                    f"percent_nuclei_{cell}_cells": float(percent),
                    f"nbr_{cell}_cells_per_mm2" : float(nbr_type_cells_per_mm2), 
                    f"nbr_nuclei_{cell}_cells" : float(nbr_nuclei_type_cells),
                    f"percent_pixels_{cell}_cells" : float(percent_pixels_type_cells),
                    "tot_cells_per_mm2" : float(tot_cells_per_mm2),
                    "tot_nuclei_wsi" : float(tot_nuclei_wsi)
                })

            results.append(CD8_dict.copy())

    # Create a dataframe with the results and reorder the columns
    results_df = pd.DataFrame(results)
    new_order_part_1 = ["Biopsi", "Group", f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2", f"nbr_nuclei_{cell}_cells", f"percent_pixels_{cell}_cells", "tot_cells_per_mm2", "tot_nuclei_wsi" ]
    new_order = new_order_part_1 + column_list
    results_df = results_df[new_order]
    
    return results_df


def before_CD8(cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_CD8: str | Path):
    
    results_df = extract_HE_and_CD8_data(column_list, file_path_HE, file_path_CD8)
    print(results_df)

    # Loop through all column names in column_list
    for i in range(len(column_list)):
        for j in range(len(HE_column_list)):
            column = column_list[i]
            column_HE = HE_column_list[j]
            
            # Calculate pearson correlation for before treatment
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b_p,p_b_p = stats.pearsonr(before_treatment[column_HE], before_treatment[column])
            r_b_s,p_b_s = stats.spearmanr(before_treatment[column_HE], before_treatment[column])

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/HE_aSMA/aSMA_{cell}/{column_HE}")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")

            plt.figure(figsize=(12, 10))
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "rebeccapurple", label =f"Pearson  P={p_b_p:.3f}, r={r_b_p:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before.png"))
            plt.show()

            plt.figure(figsize=(12, 10))
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "rebeccapurple", label =f"Spearman  P={p_b_s:.3f}, r={r_b_s:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_spearman_{column}_before.png"))
            plt.show()

def after_CD8( cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_CD8: str | Path ):
    
    results_df = extract_HE_and_CD8_data(column_list, file_path_HE, file_path_CD8)
    print(results_df)

    # Loop through all column names in column_list
    for i in range(len(column_list)):
        for j in range(len(HE_column_list)):
            column = column_list[i]
            column_HE = HE_column_list[j]
            
            # Calculate pearson correlation for after treatment
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a_p,p_a_p = stats.pearsonr(after_treatment[column_HE], after_treatment[column])
            r_a_s,p_a_s = stats.spearmanr(after_treatment[column_HE], after_treatment[column])

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/HE_aSMA/aSMA_{cell}/{column_HE}")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")

            plt.figure(figsize=(12, 10))
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "rebeccapurple", label =f"Pearson  P={p_a_p:.3f}, r={r_a_p:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_pearson_{column}_after.png"))
            plt.show()

            plt.figure(figsize=(12, 10))
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "rebeccapurple", label =f"Spearman  P={p_a_s:.3f}, r={r_a_s:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_spearman_{column}_after.png"))
            plt.show()

def before_and_after_CD8(cell: str, y_label_list: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_CD8: str | Path):
    
    results_df = extract_HE_and_CD8_data(column_list, file_path_HE, file_path_CD8)
    print(results_df)
    
    if cell=="connective":
        text_box_coord_a = [[[0,81],[5000,76]], [[60,400],[4000,7000]]] #connective
        text_box_coord_b = [[[0,69],[5000,69]], [[60,400],[4000,6200]]] #connective
        # Loop through all column names in column_list
        for i in range(len(column_list)):
            for j in range(len(HE_column_list)):
                column = column_list[i]
                column_HE = HE_column_list[j]
                
                # Calculate pearson correlation for before treatment
                before_treatment = results_df[results_df["Group"] == "before"]
                r_b_p,p_b_p = stats.pearsonr(before_treatment[column_HE], before_treatment[column])

                # Calculate pearson correlation for after treatment
                after_treatment = results_df[results_df["Group"] == "after"]
                r_a_p,p_a_p = stats.pearsonr(after_treatment[column_HE], after_treatment[column])

                
                # Path to save plots
                save_path = Path(f"/media/jenny/Expansion/MM_HE_results/correlation/pearson/HE_aSMA/aSMA_connective/all_biopsies/{column_HE}")
                if not save_path.exists():
                    save_path.mkdir(parents=True)
                    print(f"Directory {save_path} was created")

                # Fit a linear regression line
                slope_b, intercept_b, r_value, p_value, std_err = stats.linregress(before_treatment[column_HE], before_treatment[column])
                # Create a smooth x-range that goes beyond the data limits
                x_fit = np.linspace(0, before_treatment[column_HE].max()+2, 200)
                y_fit = slope_b * x_fit + intercept_b
                
                plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
                plt.figure(figsize=[12,10])
                plt.scatter(before_treatment[column_HE], before_treatment[column], color = "darkorange")
                plt.xlabel(x_label_list[j],fontsize=24, labelpad=12)
                plt.ylabel(y_label_list[i], fontsize=24, labelpad=12)
                # Tick labels
                plt.xticks(fontsize=24)
                plt.yticks(fontsize=24)
                plt.text(text_box_coord_b[i][j][0],text_box_coord_b[i][j][1], f"P = {p_b_p:.1e}, r = {r_b_p:.2f}  Before", fontsize=24, color="darkorange", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            
                # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
                # plt.legend()
                plt.plot(x_fit, y_fit , color = "darkorange")
                
                # Add labels with Biopsi IDs
                for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                    plt.text(x, y, k, fontsize=8)
                # plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before.png"))
                # plt.show()

                # Fit a linear regression line
                slope_a, intercept_a, r_value, p_value, std_err = stats.linregress(after_treatment[column_HE], after_treatment[column])
                # Create a smooth x-range that goes beyond the data limits
                x_fit = np.linspace(0, after_treatment[column_HE].max()+2, 200)
                y_fit = slope_a * x_fit + intercept_a
                
                plt.scatter(after_treatment[column_HE], after_treatment[column], color = "maroon")
                plt.xlabel(x_label_list[j], fontsize=24, labelpad=12)
                plt.ylabel(y_label_list[i], fontsize=24, labelpad=12)
                # Tick labels
                plt.xticks(fontsize=24)
                plt.yticks(fontsize=24)
                plt.text(text_box_coord_a[i][j][0],text_box_coord_a[i][j][1], f"P = {p_a_p:.1e}, r = {r_a_p:.2f}  During", fontsize=24, color="maroon", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            
                # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
                # plt.legend()
                plt.plot(x_fit, y_fit , color = "maroon")
        
                plt.title(f"Validation by IHC", fontsize=24)
                plt.tight_layout()
                # # Add labels with Biopsi IDs
                # for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
                #     plt.text(x, y, k, fontsize=8)
                plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before_and_after.png"))
                plt.show()
    
    elif cell=="lymphocyte":
        # text_box_coord_a = [[700,4500],[0,94000],[11,70],[2.5,36]] #lymphocyte
        # text_box_coord_b = [[700,4000],[0,84000],[11,62],[2.5,30]] #lymphocyte
        # text_box_coord_a = [[6000,2400]] #cells per mm2
        # text_box_coord_b = [[6000,1400]] #cells per mm2
        text_box_coord_a = [[10000,350000]] #total number of cells
        text_box_coord_b = [[10000,300000]] #total number of cells
        
        # Loop through all column names in column_list
        for i in range(len(column_list)):
            column = column_list[i]
            column_HE = HE_column_list[i]
            
            # Calculate pearson correlation for before treatment
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b_p,p_b_p = stats.pearsonr(before_treatment[column_HE], before_treatment[column])

            # Calculate pearson correlation for after treatment
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a_p,p_a_p = stats.pearsonr(after_treatment[column_HE], after_treatment[column])

            
            # Path to save plots
            # save_path = Path(f"/media/jenny/Expansion/MM_HE_results/correlation/pearson/HE_CD8/all_biopsies/{column_HE}/")
            save_path = Path(f"/media/jenny/Expansion/MM_HE_results/correlation/pearson/HE_CD8/excluding_discarded_biopsies/{column_HE}/")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")

            # Fit a linear regression line
            slope_b, intercept_b, r_value, p_value, std_err = stats.linregress(before_treatment[column_HE], before_treatment[column])
            # Create a smooth x-range that goes beyond the data limits
            x_fit = np.linspace(0, before_treatment[column_HE].max()+2, 200)
            y_fit = slope_b * x_fit + intercept_b
            
            plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
            plt.figure(figsize=[14,10])
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "darkorange")
            # start_value = 300
            # plt.ylim(bottom = start_value)
            plt.xlim(right = 350000)
            
            plt.xlabel(x_label_list[i], fontsize=24, labelpad=12)
            plt.ylabel(y_label_list[i], fontsize=24, labelpad=12)
            # Tick labels
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.text(text_box_coord_b[i][0],text_box_coord_b[i][1], f"P = {p_b_p:.1e}, r = {r_b_p:.2f}  Before", fontsize=24, color="darkorange", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
        
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            # plt.legend()
            plt.plot(x_fit, y_fit , color = "darkorange")
            
            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            # plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before.png"))
            # plt.show()

            # Fit a linear regression line
            slope_a, intercept_a, r_value, p_value, std_err = stats.linregress(after_treatment[column_HE], after_treatment[column])
            # Create a smooth x-range that goes beyond the data limits
            x_fit = np.linspace(0, after_treatment[column_HE].max()+2, 200)
            y_fit = slope_a * x_fit + intercept_a
            
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "maroon")
            plt.xlabel(x_label_list[i], fontsize=24, labelpad=12)
            plt.ylabel(y_label_list[i], fontsize=24, labelpad=12)
            # Tick labels
            plt.xticks(fontsize=24)
            plt.yticks(fontsize=24)
            plt.text(text_box_coord_a[i][0],text_box_coord_a[i][1], f"P = {p_a_p:.1e}, r = {r_a_p:.2f}  During", fontsize=24, color="maroon", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            # plt.legend()
            plt.plot(x_fit, y_fit , color = "maroon")
    
            plt.title(f"Validation by IHC", fontsize=26)
            plt.tight_layout()
            # # Add labels with Biopsi IDs
            # for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
            #     plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before_and_after.png"))
            plt.show()

def extract_gene_data_spearman_and_pearson(cell: str ,column_list: list, file_path_HE: str | Path, file_path_genes: str | Path):
    
    # Load data files as a Dataframe
    df_HE = pd.read_csv(file_path_HE)
    df_GT = pd.read_csv(file_path_genes)

    # Regex patterns for the Biopsies
    pattern_before = re.compile(r"\.[A-Z]$")
    pattern_after = re.compile(r"\.\d+[A-Z]$")

    results = []
    # Loop through all biopsies in the genes Dataframe
    for gathered_biopsies in df_GT["Biopsi"]:
        genes_dict = {}
        for column in column_list: # Loop through columns of interest to extract the data
                    col = df_GT.loc[df_GT["Biopsi"]==gathered_biopsies, column].values
                    genes_dict[column] = float(col)
        parts = gathered_biopsies.split(".")
        base = parts[0]
        subs = parts[1:]
        separate_ids = [f"{base}.{s}" for s in subs]
        
        # Group the biopsies into before and after treatment
        for biopsi in separate_ids:
            if pattern_before.search(biopsi):
                group = "before"
            elif pattern_after.search(biopsi):
                group = "after"
            else:
                None
            
            # Only include biopsies found in both input files in the results
            if biopsi in df_HE["Biopsi"].values:
                percent = df_HE.loc[df_HE["Biopsi"] == biopsi, f"percent_nuclei_{cell}_cells"].values
                print(percent)
                nbr_type_cells_per_mm2 = df_HE.loc[df_HE["Biopsi"] == biopsi, f"nbr_{cell}_cells_per_mm2"].values
                tot_cells_per_mm2 = df_HE.loc[df_HE["Biopsi"] == biopsi, "tot_cells_per_mm2"].values
                genes_dict.update({
                        "Biopsi": biopsi,
                        "Group": group,
                        f"percent_nuclei_{cell}_cells": float(percent),
                        f"nbr_{cell}_cells_per_mm2" : float(nbr_type_cells_per_mm2), 
                        "tot_cells_per_mm2" : float(tot_cells_per_mm2),
                    })
    
                results.append(genes_dict.copy())
    
    # Create a dataframe with the results and reorder the columns
    results_df = pd.DataFrame(results)
    new_order_part_1 = ["Biopsi", "Group", f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2", "tot_cells_per_mm2" ]
    new_order = new_order_part_1 + column_list
    results_df = results_df[new_order]
    
    return results_df

def extract_gene_data_spearman_and_pearson_gt_data(cell: str ,column_list: list, file_path_HE: str | Path, file_path_GT: str | Path):
    
    # Load data files as a Dataframe
    df_HE = pd.read_csv(file_path_HE)
    df_GT = pd.read_csv(file_path_GT)

    # Regex patterns for the Biopsies
    pattern_before = re.compile(r"\.[A-Z]$")
    pattern_after = re.compile(r"\.\d+[A-Z]$")

    results = []
    # Loop through all biopsies in the genes Dataframe
    for gathered_biopsies in df_GT["Biopsi"]:
        genes_dict = {}
        for column in column_list: # Loop through columns of interest to extract the data
                    col = df_GT.loc[df_GT["Biopsi"]==gathered_biopsies, column].values
                    genes_dict[column] = float(col)
        parts = gathered_biopsies.split(".")
        base = parts[0]
        subs = parts[1:]
        separate_ids = [f"{base}.{s}" for s in subs]
        
        # Group the biopsies into before and after treatment
        for biopsi in separate_ids:
            if pattern_before.search(biopsi):
                group = "before"
            elif pattern_after.search(biopsi):
                group = "after"
            else:
                None
            
            # Only include biopsies found in both input files in the results
            if biopsi in df_HE["Biopsi"].values:
                percent = df_HE.loc[df_HE["Biopsi"] == biopsi, f"percent_nuclei_{cell}_cells"].values
                print(percent)
                genes_dict.update({
                        "Biopsi": biopsi,
                        "Group": group,
                        f"percent_nuclei_{cell}_cells": float(percent),
                    })
    
                results.append(genes_dict.copy())
    
    # Create a dataframe with the results and reorder the columns
    results_df = pd.DataFrame(results)
    new_order_part_1 = ["Biopsi", "Group", f"percent_nuclei_{cell}_cells"]
    new_order = new_order_part_1 + column_list
    results_df = results_df[new_order]
    
    return results_df

def all_immune_cells_pearson_spearman(patients: int, column_list: list, file_path_genes: str | Path):
    
    # List to store dataframes from all immune cells
    dataframes = []

    stages = ["before", "after"]
    tests = ["Pearson", "Spearman"]
    x_labels = ["% immune cells", "number of cells per mm$^2$"]
    cells = ["plasma", "eosinophil", "lymphocyte", "neutrophil"]
    
    # Extract results for all immune cells
    for cell in cells:

        file_path_HE = f"/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_{cell}_HE.csv"
        results_df = extract_gene_data_spearman_and_pearson(cell, column_list, file_path_HE, file_path_genes)
        dataframes.append(results_df)
    
    headers = [" ","Test", "Measure", column_list[0]]
    test_names = []
    values = []
    signs = []
    before_after = []

    # Calculate results for all combinations of stages, tests and x_labels
    for stage, test, x_label in product(stages, tests, x_labels):
        before_after.append(stage)
        if x_label == "% immune cells":
            col_nr = 2
            name = "percentage"
            sign = "%"
            signs.append(sign)
        elif x_label == "number of cells per mm$^2$":
            col_nr = 3
            name = "nbr_per_mm2"
            sign = "mm$^2$"
            signs.append(sign)
        else:
            raise ValueError
        
        # Extract results from relevant column for before and after treatment groups
        selected_stage = [dataframe[dataframe["Group"] == stage] for dataframe in dataframes] 
        cols = [df.iloc[:, int(col_nr)] for df in selected_stage]
        
        # Extract gene information
        genes_col = selected_stage[0].iloc[:,5]
        # Add columns from all immune cells element-wise
        sum_cols = sum(cols)
        
        if test == "Pearson":
            r,p = stats.pearsonr(sum_cols, genes_col)
            test_names.append(test)
            values.append([r,p])

        elif test == "Spearman":
            r,p = stats.spearmanr(sum_cols, genes_col)
            test_names.append(test)
            values.append([r,p])
        else:
            raise ValueError

        # Path to save plots
        save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/{test}/genes_immune_cells/")
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")

        plt.figure(figsize=(12, 10))
        plt.scatter(sum_cols, genes_col, color = "rebeccapurple", label =f"{test}  P={p:.3f}, r={r:.2f}")
        plt.xlabel(x_label)
        plt.ylabel("ESTIMATE_ImmuneScore")
        plt.legend()
        
        # Add labels with Biopsi IDs
        for k, x, y in zip(selected_stage[0]["Biopsi"],sum_cols, genes_col):
            plt.text(x, y, k, fontsize=8)
        plt.savefig(save_path.joinpath(f"corr_{test}_{name}_immune_cells_{stage}.png"))
        plt.close()
    
    # Build full table rows
    table_data = []
    row_colors = []
    row_colors.append("tan") # header color

    for moment, name, sign, val in zip(before_after, test_names, signs, values):
        r, p = val
        row = [moment, name, sign, f"P = {p:.3f}, r = {r:.2f}"]

        # Row color according to significance 
        if p < 0.05:
            row_colors.append("lightgreen") 
        else:
            row_colors.append("beige")

        table_data.append(row)
    
    # Path to save tables
    save_path_tables = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/tables/genes_immune_cells/")
    if not save_path_tables.exists():
        save_path_tables.mkdir(parents=True)
        print(f"Directory {save_path_tables} was created")
    
    # Make the figure
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis("off")

    # Create the table
    table = ax.table(
        cellText=table_data,
        colLabels=headers,
        loc="center",
        cellLoc="center"
    )
    
    # Apply colors based on significance
    for i, color in enumerate(row_colors, start=0):  # include header row
        for j in range(len(table_data[0])):
            table[(i, j)].set_facecolor(color)
    
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 2)
    plt.title(f"{patients} patients")
    plt.savefig(save_path_tables.joinpath(f"immune_cells_patients_{patients}.png"))
    plt.show()
def all_immune_cells_pearson_spearman_plot_together(patients: int, column_list: list, file_path_genes: str | Path):
    
    # List to store dataframes from all immune cells
    dataframes = []

    stages = ["before", "after"]
    # tests = ["Pearson", "Spearman"]
    tests = ["Pearson"]
    x_labels = ["% Immune cells (Augmented HoVer-Net)", "Immune cells per mm$^2$ (Augmented HoVer-Net)"]
    cells = ["plasma", "eosinophil", "lymphocyte", "neutrophil"]
    
    # Extract results for all immune cells
    for cell in cells:
        file_path_HE = f"/media/jenny/Expansion1/MM_HE_results/all_statistics/all_statistics_results_{cell}_HE.csv"
        results_df = extract_gene_data_spearman_and_pearson(cell, column_list, file_path_HE, file_path_genes)
        dataframes.append(results_df)

    # Calculate results for all combinations of stages, tests and x_labels
    for test, x_label in product(tests, x_labels):

        if x_label == "% Immune cells (Augmented HoVer-Net)":
            col_nr = 2
            name = "percentage"
            sign = "%"
            text_box_coord_b = [14,2300]
            text_box_coord_a = [14,2500]
        elif x_label == "Immune cells per mm$^2$ (Augmented HoVer-Net)":
            col_nr = 3
            name = "nbr_per_mm2"
            sign = "mm$^2$"
            text_box_coord_b = [780,260]
            text_box_coord_a = [780,500]
        else:
            raise ValueError
        
        # Extract results from relevant column for before and after treatment groups
        before_stage = [dataframe[dataframe["Group"] == stages[0]] for dataframe in dataframes] 
        after_stage = [dataframe[dataframe["Group"] == stages[1]] for dataframe in dataframes]
        
        cols_b = [df.iloc[:, int(col_nr)] for df in before_stage]
        cols_a = [df.iloc[:, int(col_nr)] for df in after_stage]

        # Extract gene information
        genes_col_b = before_stage[0].iloc[:,5]
        genes_col_a = after_stage[1].iloc[:,5]
        print("------------------")
        print(len(genes_col_b))
        print(len(genes_col_a))
        # Add columns from all immune cells element-wise
        sum_cols_b = sum(cols_b)
        print(len(sum_cols_b))
        sum_cols_a = sum(cols_a)
        print(len(sum_cols_a))
        
        if test == "Pearson":
            r_b,p_b = stats.pearsonr(sum_cols_b, genes_col_b)
            r_a,p_a = stats.pearsonr(sum_cols_a, genes_col_a)

        elif test == "Spearman":
            r_b,p_b = stats.spearmanr(sum_cols_b, genes_col_b)
            r_a,p_a = stats.spearmanr(sum_cols_a, genes_col_a)
        else:
            raise ValueError

        # Path to save plots
        save_path = Path(f"/media/jenny/Expansion1/MM_HE_results/correlation/{test}/genes_immune_cells/")
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")
        
        
        # # Add labels with Biopsi IDs
        # for k, x, y in zip(selected_stage[0]["Biopsi"],sum_cols, genes_col):
        #     plt.text(x, y, k, fontsize=8)
        # plt.savefig(save_path.joinpath(f"corr_{test}_{name}_immune_cells_{stage}.png"))
        # plt.close()

        # Fit a linear regression line
        slope_a, intercept_a, r_value, p_value, std_err = stats.linregress(sum_cols_a, genes_col_a)
        # Create a smooth x-range that goes beyond the data limits
        x_fit = np.linspace(0, sum_cols_a.max()+2, 200)
        y_fit = slope_a * x_fit + intercept_a
        
        plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
        plt.figure(figsize=[12,10])
        plt.scatter(sum_cols_a, genes_col_a, color = "maroon")
        plt.xlabel(x_label, fontsize=26, labelpad=12)
        plt.ylabel("Gene-based estimate Immune score (IHC)", fontsize=26, labelpad=12)
        plt.text(text_box_coord_a[0],text_box_coord_a[1], f"P = {p_a:.3f}, r = {r_a:.2f}  During", fontsize=26, color="maroon", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
    
        # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
        # plt.legend()
        plt.plot(x_fit, y_fit , color = "maroon")

        # Fit a linear regression line
        slope_b, intercept_b, r_value, p_value, std_err = stats.linregress(sum_cols_b, genes_col_b)
        # Create a smooth x-range that goes beyond the data limits
        x_fit = np.linspace(0, sum_cols_b.max()+2, 200)
        y_fit = slope_b * x_fit + intercept_b
        
        plt.scatter(sum_cols_b, genes_col_b, color = "darkorange")
        # plt.xlabel(x_label, fontsize=21, labelpad=12)
        # plt.ylabel("Gene-based estimate Immune score", fontsize=21, labelpad=12)
        # Tick labels
        plt.xticks(fontsize=26)
        plt.yticks(fontsize=26)
        plt.text(text_box_coord_b[0],text_box_coord_b[1], f"P = {p_b:.3f}, r = {r_b:.2f}  Before", fontsize=26, color="darkorange", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
    
        # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
        # plt.legend()
        plt.plot(x_fit, y_fit , color = "darkorange")
        plt.title(f"Validation by gene-based estimate", fontsize=26)
        plt.tight_layout()
        plt.savefig(save_path.joinpath(f"corr_{test}_{name}_immune_cells.png"))
        plt.show()


def pearson_test(patients: int, cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_genes: str | Path):           
    
    results_df = extract_gene_data_spearman_and_pearson(cell, column_list, file_path_HE, file_path_genes)
    row_colors = []
    row_colors.append("tan") # header color

    # Prepare dict: shared columns + one list per data column
    column_dict = {"Time": [], "Test": [], "Measure": []}
    color_dict = {}
    for col in column_list:
        column_dict[col] = []
        color_dict[col] = []

    if cell == "epithelial":
        text_box_coord_b = [[[48, 1250]],[[3100, 1250]]] #epithelial
        text_box_coord_a = [[[48, 1550]],[[3100, 1550]]] #epithelial
    elif cell == "lymphocyte":
        text_box_coord_b = [[[14, 6.8],[14, 2180],[14, 5.25], [14, 7.6], [14, 7.7], [14, 4.2]],[[14, 6.8],[14, 2180],[14, 5.25], [14, 7.6], [14, 7.7], [750, 2.32]]] #lymphocyte
        text_box_coord_a = [[[14, 7.2],[14, 2500],[14, 5.55], [14, 8.1], [14, 8.2], [14, 4.55]], [[14, 7.2],[14, 2500],[14, 5.55], [14, 8.1], [14, 8.2], [750, 2.55]]] #lymphocyte
    else:
        print(f"No textbox coordinates created for {cell}")

    # For each HE column create exactly ONE row (shared entries) and append each column's value
    for i, column_HE in enumerate(HE_column_list):
        # append shared fields once per HE column
        column_dict["Time"].append("before")
        column_dict["Test"].append("Pearson")

        if column_HE == f"percent_nuclei_{cell}_cells":
            column_dict["Measure"].append("%")
            x_label = x_label_list[0]
        elif column_HE == f"nbr_{cell}_cells_per_mm2":
            column_dict["Measure"].append("mm$^2$")
            x_label = x_label_list[1]
        else:
            raise ValueError(f"Unexpected HE column: {column_HE}")

        # For each requested output column, compute the r/p and append to that column's list
        for j, column in enumerate(column_list):
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b, p_b = stats.pearsonr(before_treatment[column_HE], before_treatment[column])
            column_dict[column].append(f"P = {p_b:.1e}, r = {r_b:.2f}")
            # Row color according to significance 
            if p_b < 0.05:
                color_dict[column].append("lightgreen") 
            else:
                color_dict[column].append("beige")

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion1/MM_HE_results/correlation/pearson/genes_{cell}/{column_HE}")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")
            
            if cell == "epithelial":
                y_label = "Gene-based estimate Stromal score (IHC)"
            elif cell == "lymphocyte":
                y_label = column
            else:
                print(f"No y-label created for {cell}")

            # Fit a linear regression line
            slope_b, intercept_b, r_value, p_value, std_err = stats.linregress(before_treatment[column_HE], before_treatment[column])
            # Create a smooth x-range that goes beyond the data limits
            x_fit = np.linspace(0, before_treatment[column_HE].max()+2, 200)
            y_fit = slope_b * x_fit + intercept_b
            plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
            plt.figure(figsize=[12,10])
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "darkorange")
            plt.text(text_box_coord_b[i][j][0],text_box_coord_b[i][j][1], f"P = {p_b:.1e}, r = {r_b:.2f}  Before", fontsize=26, color="darkorange", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.plot(x_fit, y_fit , color = "darkorange")
            # plt.xlabel(x_label)
            # plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            # plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            # plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before_poster_no_labels.png"))
            # # plt.close()
            # plt.show()

        # column_dict["Time"].append("after")
        # column_dict["Test"].append("Pearson")

        # if column_HE == f"percent_nuclei_{cell}_cells":
        #     column_dict["Measure"].append("%")
        #     x_label = x_label_list[0]
        # elif column_HE == f"nbr_{cell}_cells_per_mm2":
        #     column_dict["Measure"].append("mm$^2$")
        #     x_label = x_label_list[1]
        # else:
        #     raise ValueError(f"Unexpected HE column: {column_HE}")

        # For each requested output column, compute the r/p and append to that column's list
        # for column in column_list:
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a, p_a = stats.pearsonr(after_treatment[column_HE], after_treatment[column])
            column_dict[column].append(f"P = {p_a:.1e}, r = {r_a:.2f}")
            # Row color according to significance 
            if p_a < 0.05:
                color_dict[column].append("lightgreen") 
            else:
                color_dict[column].append("beige")
            
            # Fit a linear regression line
            slope_a, intercept_a, r_value, p_value, std_err = stats.linregress(after_treatment[column_HE], after_treatment[column])
            # Create a smooth x-range that goes beyond the data limits
            x_fit = np.linspace(0, after_treatment[column_HE].max()+2, 200)
            y_fit = slope_a * x_fit + intercept_a
            
            # plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
            # plt.figure(figsize=[12,10])
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "maroon")
            plt.text(text_box_coord_a[i][j][0],text_box_coord_a[i][j][1], f"P = {p_a:.1e}, r = {r_a:.2f}  During", fontsize=26, color="maroon", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.xlabel(x_label, fontsize=26, labelpad=12)
            plt.ylabel(y_label, fontsize=26, labelpad=12)
            # Tick labels
            plt.xticks(fontsize=26)
            plt.yticks(fontsize=26)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            # plt.legend()
            plt.plot(x_fit, y_fit , color = "maroon")
            # # Add labels with Biopsi IDs
            # for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
            #     plt.text(x, y, k, fontsize=8)
            plt.title(f"Validation by gene-based estimate", fontsize=26)
            plt.tight_layout()   
            plt.savefig(save_path.joinpath(f"corr_pearson_{column}.png"))
            plt.show()
            # plt.close()

    
    # df = pd.DataFrame(column_dict)
    # print(df)

    # # Build a cell_colors 2D list for the table
    # cell_colors = []

    # for i in range(len(df)):
    #     row_colors = []
    #     for col in df.columns:
    #         # Only apply colors to columns that exist in color_dict, otherwise white
    #         if col in color_dict:
    #             row_colors.append(color_dict[col][i])
    #         else:
    #             row_colors.append("beige")  # shared columns like Time/Test/Sign
    #     cell_colors.append(row_colors)
    
    # # Path to save tables
    # save_path_tables = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/tables/genes_{cell}/")
    # if not save_path_tables.exists():
    #     save_path_tables.mkdir(parents=True)
    #     print(f"Directory {save_path_tables} was created")
    
    # # Create figure
    # fig, ax = plt.subplots(figsize=(18, 4))  # adjust size to fit table
    # ax.axis("off")  # hide axes

    # # Create table
    # tbl = ax.table(cellText=df.values, colLabels=df.columns, cellColours=cell_colors, loc='center', cellLoc='center')
    # # Set color of header
    # for j in range(len(df.columns)):
    #     tbl[(0, j)].set_facecolor("tan")  # header row is row 0
    # tbl.auto_set_font_size(False)
    # tbl.set_fontsize(12)
    # tbl.auto_set_column_width(col=list(range(len(df.columns))))
    # tbl.scale(1.2, 2)
    # plt.title(f"{patients} patients")
    # plt.savefig(save_path_tables.joinpath(f"{cell}_pearson_patients_{patients}.png"))
    # plt.plot()
    
def pearson_test_gt_data(patients: int, cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_GT: str | Path):           
    
    results_df = extract_gene_data_spearman_and_pearson_gt_data(cell, column_list, file_path_HE, file_path_GT)
    row_colors = []
    row_colors.append("tan") # header color

    # Prepare dict: shared columns + one list per data column
    column_dict = {"Time": [], "Test": [], "Measure": []}
    color_dict = {}
    for col in column_list:
        column_dict[col] = []
        color_dict[col] = []

    if cell == "epithelial":
        text_box_coord_b = [[[0, 95]]] #epithelial
        text_box_coord_a = [[[0, 83]]] #epithelial
    elif cell == "lymphocyte":
        text_box_coord_b = [[[14, 6.8],[14, 2180],[14, 5.25], [14, 7.6], [14, 7.7], [14, 4.2]],[[14, 6.8],[14, 2180],[14, 5.25], [14, 7.6], [14, 7.7], [750, 2.32]]] #lymphocyte
        text_box_coord_a = [[[14, 7.2],[14, 2500],[14, 5.55], [14, 8.1], [14, 8.2], [14, 4.55]], [[14, 7.2],[14, 2500],[14, 5.55], [14, 8.1], [14, 8.2], [750, 2.55]]] #lymphocyte
    else:
        print(f"No textbox coordinates created for {cell}")

    # For each HE column create exactly ONE row (shared entries) and append each column's value
    for i, column_HE in enumerate(HE_column_list):
        # append shared fields once per HE column
        column_dict["Time"].append("before")
        column_dict["Test"].append("Pearson")

        if column_HE == f"percent_nuclei_{cell}_cells":
            column_dict["Measure"].append("%")
            x_label = x_label_list[0]
        else:
            raise ValueError(f"Unexpected HE column: {column_HE}")

        # For each requested output column, compute the r/p and append to that column's list
        for j, column in enumerate(column_list):
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b, p_b = stats.pearsonr(before_treatment[column_HE], before_treatment[column])
            column_dict[column].append(f"P = {p_b:.1e}, r = {r_b:.2f}")
            # Row color according to significance 
            if p_b < 0.05:
                color_dict[column].append("lightgreen") 
            else:
                color_dict[column].append("beige")

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/correlation/pearson/{cell}/GT_and_hovernet/")
        
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")
            
            if cell == "epithelial":
                y_label = "Ground truth tumor percentage"
            elif cell == "lymphocyte":
                y_label = column
            else:
                print(f"No y-label created for {cell}")

            # Fit a linear regression line
            slope_b, intercept_b, r_value, p_value, std_err = stats.linregress(before_treatment[column_HE], before_treatment[column])
            # Create a smooth x-range that goes beyond the data limits
            x_fit = np.linspace(0, before_treatment[column_HE].max()+2, 200)
            y_fit = slope_b * x_fit + intercept_b
            plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
            plt.figure(figsize=[12,10])
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "darkorange")
            plt.text(text_box_coord_b[i][j][0],text_box_coord_b[i][j][1], f"P = {p_b:.1e}, r = {r_b:.2f}  Before", fontsize=26, color="darkorange", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.plot(x_fit, y_fit , color = "darkorange")
            # plt.xlabel(x_label)
            # plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            # plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            # plt.savefig(save_path.joinpath(f"corr_pearson_{column}_before_poster_no_labels.png"))
            # # plt.close()
            # plt.show()

        # column_dict["Time"].append("after")
        # column_dict["Test"].append("Pearson")

        # if column_HE == f"percent_nuclei_{cell}_cells":
        #     column_dict["Measure"].append("%")
        #     x_label = x_label_list[0]
        # elif column_HE == f"nbr_{cell}_cells_per_mm2":
        #     column_dict["Measure"].append("mm$^2$")
        #     x_label = x_label_list[1]
        # else:
        #     raise ValueError(f"Unexpected HE column: {column_HE}")

        # For each requested output column, compute the r/p and append to that column's list
        # for column in column_list:
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a, p_a = stats.pearsonr(after_treatment[column_HE], after_treatment[column])
            column_dict[column].append(f"P = {p_a:.1e}, r = {r_a:.2f}")
            # Row color according to significance 
            if p_a < 0.05:
                color_dict[column].append("lightgreen") 
            else:
                color_dict[column].append("beige")
            
            # Fit a linear regression line
            slope_a, intercept_a, r_value, p_value, std_err = stats.linregress(after_treatment[column_HE], after_treatment[column])
            # Create a smooth x-range that goes beyond the data limits
            x_fit = np.linspace(0, after_treatment[column_HE].max()+2, 200)
            y_fit = slope_a * x_fit + intercept_a
            
            # plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
            # plt.figure(figsize=[12,10])
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "maroon")
            plt.text(text_box_coord_a[i][j][0],text_box_coord_a[i][j][1], f"P = {p_a:.1e}, r = {r_a:.2f}  During", fontsize=26, color="maroon", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.xlabel(x_label, fontsize=26, labelpad=12)
            plt.ylabel(y_label, fontsize=26, labelpad=12)
            # Tick labels
            plt.xticks(fontsize=26)
            plt.yticks(fontsize=26)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            # plt.legend()
            plt.plot(x_fit, y_fit , color = "maroon")
            # # Add labels with Biopsi IDs
            # for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
            #     plt.text(x, y, k, fontsize=8)
            plt.title(f"Validation by ground truth", fontsize=26)
            plt.tight_layout()   
            plt.savefig(save_path.joinpath(f"corr_pearson_{column}.png"))
            plt.show()
            # plt.close()

    
    # df = pd.DataFrame(column_dict)
    # print(df)

    # # Build a cell_colors 2D list for the table
    # cell_colors = []

    # for i in range(len(df)):
    #     row_colors = []
    #     for col in df.columns:
    #         # Only apply colors to columns that exist in color_dict, otherwise white
    #         if col in color_dict:
    #             row_colors.append(color_dict[col][i])
    #         else:
    #             row_colors.append("beige")  # shared columns like Time/Test/Sign
    #     cell_colors.append(row_colors)
    
    # # Path to save tables
    # save_path_tables = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/tables/genes_{cell}/")
    # if not save_path_tables.exists():
    #     save_path_tables.mkdir(parents=True)
    #     print(f"Directory {save_path_tables} was created")
    
    # # Create figure
    # fig, ax = plt.subplots(figsize=(18, 4))  # adjust size to fit table
    # ax.axis("off")  # hide axes

    # # Create table
    # tbl = ax.table(cellText=df.values, colLabels=df.columns, cellColours=cell_colors, loc='center', cellLoc='center')
    # # Set color of header
    # for j in range(len(df.columns)):
    #     tbl[(0, j)].set_facecolor("tan")  # header row is row 0
    # tbl.auto_set_font_size(False)
    # tbl.set_fontsize(12)
    # tbl.auto_set_column_width(col=list(range(len(df.columns))))
    # tbl.scale(1.2, 2)
    # plt.title(f"{patients} patients")
    # plt.savefig(save_path_tables.joinpath(f"{cell}_pearson_patients_{patients}.png"))
    # plt.plot()
def spearman_test(patients: int, cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_genes: str | Path):           
    
    results_df = extract_gene_data_spearman_and_pearson(cell, column_list, file_path_HE, file_path_genes)
    row_colors = []
    row_colors.append("tan") # header color

    # Prepare dict: shared columns + one list per data column
    column_dict = {"Time": [], "Test": [], "Measure": []}
    color_dict = {}
    for col in column_list:
        column_dict[col] = []
        color_dict[col] = []
    
    # For each HE column create exactly ONE row (shared entries) and append each column's value
    for column_HE in HE_column_list:
        # append shared fields once per HE column
        column_dict["Time"].append("before")
        column_dict["Test"].append("Spearman")

        if column_HE == f"percent_nuclei_{cell}_cells":
            column_dict["Measure"].append("%")
            x_label = x_label_list[0]
        elif column_HE == f"nbr_{cell}_cells_per_mm2":
            column_dict["Measure"].append("mm$^2$")
            x_label = x_label_list[1]
        else:
            raise ValueError(f"Unexpected HE column: {column_HE}")

        # For each requested output column, compute the r/p and append to that column's list
        for column in column_list:
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b, p_b = stats.spearmanr(before_treatment[column_HE], before_treatment[column])
            column_dict[column].append(f"P = {p_b:.3f}, r = {r_b:.2f}")
            # Row color according to significance 
            if p_b < 0.05:
                color_dict[column].append("lightgreen") 
            else:
                color_dict[column].append("beige")

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/spearman/genes_{cell}/{column_HE}")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")

            plt.figure(figsize=(12, 10))
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "darkorange", label =f"Spearman  P={p_b:.3f}, r={r_b:.2f}")
            plt.xlabel(x_label)
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            # plt.savefig(save_path.joinpath(f"corr_spearman_{column}_before.png"))
            plt.close()

        column_dict["Time"].append("after")
        column_dict["Test"].append("Spearman")

        if column_HE == f"percent_nuclei_{cell}_cells":
            column_dict["Measure"].append("%")
            x_label = x_label_list[0]
        elif column_HE == f"nbr_{cell}_cells_per_mm2":
            column_dict["Measure"].append("mm$^2$")
            x_label = x_label_list[1]
        else:
            raise ValueError(f"Unexpected HE column: {column_HE}")

        # For each requested output column, compute the r/p and append to that column's list
        for column in column_list:
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a, p_a = stats.spearmanr(after_treatment[column_HE], after_treatment[column])
            column_dict[column].append(f"P = {p_a:.3f}, r = {r_a:.2f}")
            # Row color according to significance 
            if p_a < 0.05:
                color_dict[column].append("lightgreen") 
            else:
                color_dict[column].append("beige")

            plt.figure(figsize=(12, 10))
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "rebeccapurple", label =f"Spearman  P={p_a:.3f}, r={r_a:.2f}")
            plt.xlabel(x_label)
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            # plt.savefig(save_path.joinpath(f"corr_spearman_{column}_after.png"))
            plt.close()

    
    # df = pd.DataFrame(column_dict)
    # print(df)

    # # Build a cell_colors 2D list for the table
    # cell_colors = []

    # for i in range(len(df)):
    #     row_colors = []
    #     for col in df.columns:
    #         # Only apply colors to columns that exist in color_dict, otherwise white
    #         if col in color_dict:
    #             row_colors.append(color_dict[col][i])
    #         else:
    #             row_colors.append("beige")  # shared columns like Time/Test/Sign
    #     cell_colors.append(row_colors)
    
    # # Path to save tables
    # save_path_tables = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/tables/genes_{cell}/")
    # if not save_path_tables.exists():
    #     save_path_tables.mkdir(parents=True)
    #     print(f"Directory {save_path_tables} was created")
    
    # # Create figure
    # fig, ax = plt.subplots(figsize=(18, 4))  # adjust size to fit table
    # ax.axis("off")  # hide axes

    # # Create table
    # tbl = ax.table(cellText=df.values, colLabels=df.columns, cellColours=cell_colors, loc='center', cellLoc='center')
    # # Set color of header
    # for j in range(len(df.columns)):
    #     tbl[(0, j)].set_facecolor("tan")  # header row is row 0
    # tbl.auto_set_font_size(False)
    # tbl.set_fontsize(12)
    # tbl.auto_set_column_width(col=list(range(len(df.columns))))
    # tbl.scale(1.2, 2)
    # plt.title(f"{patients} patients")
    # plt.savefig(save_path_tables.joinpath(f"{cell}_spearman_patients_{patients}_poster.png"))
    # plt.plot()

def extract_data_for_wilcoxon_and_t_test(cell: str, file_path_HE: str | Path):

    df = pd.read_csv(file_path_HE)

    # Regex patterns for the Biopsies
    pattern_before = re.compile(r"\.[A-Z]$")
    pattern_after = re.compile(r"\.\d+[A-Z]$")

    # Extract base ID (NN041, NN201, ...) from Biopsi column and add to new "base" column
    df["base"] = df["Biopsi"].str.extract(r"^(.*?)(?=\.\d*[A-Z]$)")

    # Classify before and after
    df["group"] = df["Biopsi"].apply(
        lambda x: "before" if pattern_before.search(x)
                else "after" if pattern_after.search(x)
                else None
    )

    results = []
    # Extract the biopsies that belong to before and after
    for base, group in df.groupby("base"):
        befores = group[group["group"] == "before"]["Biopsi"].unique()
        afters  = group[group["group"] == "after"]["Biopsi"].unique()
        for b in befores:
            for a in afters:
                # Collect all replicates for this before and after pair
                percent_before = df.loc[df["Biopsi"] == b, f"percent_nuclei_{cell}_cells"].values
                percent_after = df.loc[df["Biopsi"] == a, f"percent_nuclei_{cell}_cells"].values
        
                nbr_mm2_before = df.loc[df["Biopsi"] == b, f"nbr_{cell}_cells_per_mm2"].values
                nbr_mm2_after = df.loc[df["Biopsi"] == a, f"nbr_{cell}_cells_per_mm2"].values
                
                tot_cells_per_mm2_before = df.loc[df["Biopsi"] == b, "tot_cells_per_mm2"].values
                tot_cells_per_mm2_after = df.loc[df["Biopsi"] == a, "tot_cells_per_mm2"].values
                
                results.append({
                    "base": base,
                    "before": b,
                    "after": a,
                    "percent_before": percent_before,
                    "percent_after": percent_after,
                    "nbr_mm2_before": nbr_mm2_before,
                    "nbr_mm2_after": nbr_mm2_after,
                    "tot_cells_per_mm2_before": tot_cells_per_mm2_before,
                    "tot_cells_per_mm2_after": tot_cells_per_mm2_after
                })
    # Create a dataframe with the results
    results_df = pd.DataFrame(results)
    return results_df


def wilcoxon_test(column_list: list, y_label_list: list, cell: str, file_path_HE: str | Path):
    
    # Load dataframe
    results_df = extract_data_for_wilcoxon_and_t_test(cell, file_path_HE)  
    
    # Convert to cells per micrometer (Log)
    results_df["nbr_mm2_before"] = results_df["nbr_mm2_before"].apply(lambda x: x[0]) 
    results_df["nbr_mm2_before"] = np.log(results_df["nbr_mm2_before"]/(1e6))
    
    results_df["nbr_mm2_after"] = results_df["nbr_mm2_after"].apply(lambda x: x[0])
    results_df["nbr_mm2_after"] = np.log(results_df["nbr_mm2_after"]/(1e6))
    
    results_df["tot_cells_per_mm2_before"] = results_df["tot_cells_per_mm2_before"].apply(lambda x: x[0])
    results_df["tot_cells_per_mm2_before"] = np.log(results_df["tot_cells_per_mm2_before"]/(1e6))
    
    results_df["tot_cells_per_mm2_after"] = results_df["tot_cells_per_mm2_after"].apply(lambda x: x[0])
    results_df["tot_cells_per_mm2_after"] = np.log(results_df["tot_cells_per_mm2_after"]/(1e6))

    if cell == "neutrophil":
        text_box_pos = [[1.22,0.6], [1.22, -11], [1.22, -4.5]]  #Neutrophil
    elif cell == "eosinophil":
        text_box_pos = [[1.22,0.05], [1.22, -12.5], [1.22, -4.5]]  #Eosinophil
    elif cell == "connective":
        text_box_pos = [[1.232,98], [1.22, -5], [1.22, -4.5]] #connective
    elif cell == "epithelial":
        text_box_pos = [[1.22,80], [1.22, -5], [1.22, 11000]] #epithel
    elif cell == "plasma":
        text_box_pos = [[1.22,5.7], [1.22, -8], [1.22, -4.5]]  #plasma
    elif cell == "lymphocyte":
        text_box_pos = [[1.22,17.5], [1.22, -7], [1.22, -4.5]]  #lymphocyte
    else:
        print("No valid cell type chosen")
    
    # Loop through all column names in column_list
    for i in range(len(column_list)):
        column = column_list[i]
        y_label = y_label_list[i]
        print(results_df[column + "_before"])
        before_treatment = results_df[column + "_before"].to_numpy(dtype=float)
        after_treatment = results_df[column + "_after"].to_numpy(dtype=float)
        # print(before_treatment)
        # print(after_treatment)
        if before_treatment.size != after_treatment.size:
            print("before_treatment and after_treatment is not of same size")
        
        # Calculate the Wilcoxon signed-rank test
        try:
            stat, p = wilcoxon(before_treatment, after_treatment)
        
        except ValueError:
            stat, p = None, None
            print("p = None due to ValueError")

        print(f"Wilcoxon signed-rank test on {y_label} gives: stat={stat},  p={p}")
        
        # X positions for the bars
        x_positions = [1, 1.2]
        
        plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
        plt.figure(figsize=[10,8])

        # # Bars
        # plt.bar(x_positions, max_vals, color=['darkorange', 'darkblue'], alpha=0.6, width=0.6)
        box = plt.boxplot(
                    [before_treatment, after_treatment],
                    positions=x_positions,
                    patch_artist=True,    # fill the box with color
                    medianprops=dict(color='black', linewidth=2),
                    # whiskerprops=dict(color='gray', linewidth=1.5),
                    capprops=dict(linewidth=0)  # removes the little line at whisker ends
                )
        
        # Box colors
        box_colors = ["darkorange", "maroon"]
        for patch, color in zip(box["boxes"], box_colors):
            patch.set_facecolor(color)
        # Overlay dot points
        # Add small random jitter so dots don’t overlap exactly
        jitter_before = np.random.normal(x_positions[0], 0.04, len(before_treatment))
        jitter_after = np.random.normal(x_positions[1], 0.04, len(after_treatment))
        x_pos_array_before = np.full(len(before_treatment), x_positions[0])
        x_pos_array_after = np.full(len(after_treatment), x_positions[1])
        
        plt.scatter(x_pos_array_before, before_treatment, alpha=0.6, color='black', s=20, zorder=3)
        plt.scatter(x_pos_array_after, after_treatment, alpha=0.6, color='black', s=20, zorder=3)
        
        # Connect each dot pair with a line
        for xb, yb, xa, ya in zip(x_pos_array_before, before_treatment, x_pos_array_after, after_treatment):
            plt.plot([xb, xa], [yb, ya], color='gray', alpha=0.5, linewidth=1, zorder=2)

        save_path = Path(f"/media/jenny/Expansion/MM_HE_results/correlation/wilcoxon/excluding_discarded_biopsies/{cell}/micrometer/")
        
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")

        # Add labels
        plt.text(text_box_pos[i][0], text_box_pos[i][1], f"P = {p:.1e}", fontsize=21, color="black", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
        plt.xticks(x_positions, ["Before", "During"], fontsize=21)
        plt.yticks(fontsize=21)
        plt.ylabel(y_label, fontsize=21, labelpad=12)
        plt.xlim(0.85, 1.35)
        plt.title(f"HoVer-Net", fontsize=21)
        plt.tight_layout()
        plt.savefig(save_path.joinpath(f"wilcoxon_{column}.png"))
        plt.show()

def t_test(column_list: list, y_label_list: list, cell: str, file_path_HE: str | Path):
    
    # Load dataframe
    results_df = extract_data_for_wilcoxon_and_t_test(cell, file_path_HE)
    
    epsilon_b = 1e-11
    epsilon_a = 1e-10
    # Convert to cells per micrometer (Log)
    results_df["nbr_mm2_before"] = results_df["nbr_mm2_before"].apply(lambda x: x[0]) 
    results_df["nbr_mm2_before"] = np.log((results_df["nbr_mm2_before"]/(1e6)) + epsilon_b)
    
    results_df["nbr_mm2_after"] = results_df["nbr_mm2_after"].apply(lambda x: x[0])
    results_df["nbr_mm2_after"] = np.log((results_df["nbr_mm2_after"]/(1e6)) + epsilon_a)
    
    results_df["tot_cells_per_mm2_before"] = results_df["tot_cells_per_mm2_before"].apply(lambda x: x[0])
    results_df["tot_cells_per_mm2_before"] = np.log(results_df["tot_cells_per_mm2_before"]/(1e6))
    
    results_df["tot_cells_per_mm2_after"] = results_df["tot_cells_per_mm2_after"].apply(lambda x: x[0])
    results_df["tot_cells_per_mm2_after"] = np.log(results_df["tot_cells_per_mm2_after"]/(1e6))

    if cell == "neutrophil":
        text_box_pos = [[1.22,0.6], [1.22, -11], [1.22, -4.5]]  #Neutrophil
    elif cell == "eosinophil":
        text_box_pos = [[1.22,0.05], [1.22, -12.5], [1.22, -4.5]]  #Eosinophil
    elif cell == "connective":
        text_box_pos = [[1.232,98], [1.22, -5], [1.22, -4.5]] #connective
    elif cell == "epithelial":
        text_box_pos = [[1.22,80], [1.22, -5], [1.22, -4.5]] #epithel
    elif cell == "plasma":
        text_box_pos = [[1.22,5.7], [1.22, -8], [1.22, -4.5]]  #plasma
    elif cell == "lymphocyte":
        text_box_pos = [[1.22,17.5], [1.22, -7], [1.22, -4.5]]  #lymphocyte
    else:
        print("No valid cell type chosen")
    
    # Loop through all column names in column_list
    for i in range(len(column_list)):
        column = column_list[i]
        y_label = y_label_list[i]
        print(results_df[column + "_before"])
        before_treatment = results_df[column + "_before"].to_numpy(dtype=float)
        after_treatment = results_df[column + "_after"].to_numpy(dtype=float)
        # print(before_treatment)
        # print(after_treatment)
        if before_treatment.size != after_treatment.size:
            print("before_treatment and after_treatment is not of same size")
        
        diff = before_treatment - after_treatment
        print("diff array:", diff)
        print("diff variance:", np.var(diff))
        # Perform t-test
        try:
            stat, p = ttest_rel(before_treatment, after_treatment)
        
        except ValueError:
            stat, p = None, None
            print("p = None due to ValueError")
       
        print(f"t-test on {y_label} gives: stat={stat},  p={p}")
        
        # X positions for the bars
        x_positions = [1, 1.2]
        
        plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
        plt.figure(figsize=[10,8])

        # # Bars
        # plt.bar(x_positions, max_vals, color=['darkorange', 'darkblue'], alpha=0.6, width=0.6)
        box = plt.boxplot(
                    [before_treatment, after_treatment],
                    positions=x_positions,
                    patch_artist=True,    # fill the box with color
                    medianprops=dict(color='black', linewidth=2),
                    # whiskerprops=dict(color='gray', linewidth=1.5),
                    capprops=dict(linewidth=0)  # removes the little line at whisker ends
                )
        
        # Box colors
        box_colors = ["darkorange", "maroon"]
        for patch, color in zip(box["boxes"], box_colors):
            patch.set_facecolor(color)
        # Overlay dot points
        # Add small random jitter so dots don’t overlap exactly
        jitter_before = np.random.normal(x_positions[0], 0.04, len(before_treatment))
        jitter_after = np.random.normal(x_positions[1], 0.04, len(after_treatment))
        x_pos_array_before = np.full(len(before_treatment), x_positions[0])
        x_pos_array_after = np.full(len(after_treatment), x_positions[1])
        
        plt.scatter(x_pos_array_before, before_treatment, alpha=0.6, color='black', s=20, zorder=3)
        plt.scatter(x_pos_array_after, after_treatment, alpha=0.6, color='black', s=20, zorder=3)
        
        # Connect each dot pair with a line
        for xb, yb, xa, ya in zip(x_pos_array_before, before_treatment, x_pos_array_after, after_treatment):
            plt.plot([xb, xa], [yb, ya], color='gray', alpha=0.5, linewidth=1, zorder=2)

        save_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/correlation/t_test/{cell}/micrometer/")
        
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")

        # Add labels
        plt.text(text_box_pos[i][0], text_box_pos[i][1], f"P = {p:.1e}", fontsize=21, color="black", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
        plt.xticks(x_positions, ["Before", "During"], fontsize=21)
        plt.yticks(fontsize=21)
        plt.ylabel(y_label, fontsize=21, labelpad=12)
        plt.xlim(0.85, 1.35)
        plt.title(f"HoVer-Net", fontsize=21)
        plt.tight_layout()
        plt.savefig(save_path.joinpath(f"t_test_{column}.png"))
        plt.show()

def extract_data_for_wilcoxon_and_t_test_gt_data(cell: str, file_path_HE: str | Path):

    df = pd.read_csv(file_path_HE)

    # Regex patterns for the Biopsies
    pattern_before = re.compile(r"\.[A-Z]$")
    pattern_after = re.compile(r"\.\d+[A-Z]$")

    # Extract base ID (NN041, NN201, ...) from Biopsi column and add to new "base" column
    df["base"] = df["Biopsi"].str.extract(r"^(.*?)(?=\.\d*[A-Z]$)")

    # Classify before and after
    df["group"] = df["Biopsi"].apply(
        lambda x: "before" if pattern_before.search(x)
                else "after" if pattern_after.search(x)
                else None
    )

    results = []
    # Extract the biopsies that belong to before and after
    for base, group in df.groupby("base"):
        befores = group[group["group"] == "before"]["Biopsi"].unique()
        afters  = group[group["group"] == "after"]["Biopsi"].unique()
        for b in befores:
            for a in afters:
                # Collect all replicates for this before and after pair
                percent_before = df.loc[df["Biopsi"] == b, f"percentage_tumor_Heide"].values
                percent_after = df.loc[df["Biopsi"] == a, f"percentage_tumor_Heide"].values
        
                
                results.append({
                    "base": base,
                    "before": b,
                    "after": a,
                    "percent_before": percent_before,
                    "percent_after": percent_after
                })
    # Create a dataframe with the results
    results_df = pd.DataFrame(results)
    return results_df


def t_test_gt_data(column_list: list, y_label_list: list, cell: str, file_path_HE: str | Path):
    
    # Load dataframe
    results_df = extract_data_for_wilcoxon_and_t_test_gt_data(cell, file_path_HE)
    
    # Loop through all column names in column_list
    for i in range(len(column_list)):
        column = column_list[i]
        y_label = y_label_list[i]
        print(results_df[column + "_before"])
        before_treatment = results_df[column + "_before"].to_numpy(dtype=float)
        after_treatment = results_df[column + "_after"].to_numpy(dtype=float)
        # print(before_treatment)
        # print(after_treatment)
        if before_treatment.size != after_treatment.size:
            print("before_treatment and after_treatment is not of same size")
        
        diff = before_treatment - after_treatment
        print("diff array:", diff)
        print("diff variance:", np.var(diff))
        # Perform t-test
        try:
            stat, p = ttest_rel(before_treatment, after_treatment)
        
        except ValueError:
            stat, p = None, None
            print("p = None due to ValueError")
       
        print(f"t-test on {y_label} gives: stat={stat},  p={p}")
        
        # X positions for the bars
        x_positions = [1, 1.2]
        
        plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
        plt.figure(figsize=[10,8])

        # # Bars
        # plt.bar(x_positions, max_vals, color=['darkorange', 'darkblue'], alpha=0.6, width=0.6)
        box = plt.boxplot(
                    [before_treatment, after_treatment],
                    positions=x_positions,
                    patch_artist=True,    # fill the box with color
                    medianprops=dict(color='black', linewidth=2),
                    # whiskerprops=dict(color='gray', linewidth=1.5),
                    capprops=dict(linewidth=0)  # removes the little line at whisker ends
                )
        
        # Box colors
        box_colors = ["darkorange", "maroon"]
        for patch, color in zip(box["boxes"], box_colors):
            patch.set_facecolor(color)
        # Overlay dot points
        # Add small random jitter so dots don’t overlap exactly
        jitter_before = np.random.normal(x_positions[0], 0.04, len(before_treatment))
        jitter_after = np.random.normal(x_positions[1], 0.04, len(after_treatment))
        x_pos_array_before = np.full(len(before_treatment), x_positions[0])
        x_pos_array_after = np.full(len(after_treatment), x_positions[1])
        
        plt.scatter(x_pos_array_before, before_treatment, alpha=0.6, color='black', s=20, zorder=3)
        plt.scatter(x_pos_array_after, after_treatment, alpha=0.6, color='black', s=20, zorder=3)
        
        # Connect each dot pair with a line
        for xb, yb, xa, ya in zip(x_pos_array_before, before_treatment, x_pos_array_after, after_treatment):
            plt.plot([xb, xa], [yb, ya], color='gray', alpha=0.5, linewidth=1, zorder=2)

        save_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/correlation/t_test/{cell}/GT_only_before_after/")
        
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")
        
        text_box_pos = [[1.22,90], [1.22, -5], [1.22, -4.5]] 

        # Add labels
        plt.text(text_box_pos[i][0], text_box_pos[i][1], f"P = {p:.1e}", fontsize=21, color="black", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
        plt.xticks(x_positions, ["Before", "During"], fontsize=21)
        plt.yticks(fontsize=21)
        plt.ylabel(y_label, fontsize=21, labelpad=12)
        plt.xlim(0.85, 1.35)
        plt.title(f"Ground truth", fontsize=21)
        plt.tight_layout()
        plt.savefig(save_path.joinpath(f"t_test_{column}.png"))
        plt.show()

def cell_ratio(cells):
    
    epsilon = 1e-10
    # Load dataframes
    file_path_HE_1 = f"/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_results_{cells[0]}_HE.csv"
    file_path_HE_2 = f"/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_results_{cells[1]}_HE.csv"
    
    results_df_1 = extract_data_for_wilcoxon_and_t_test(cells[0], file_path_HE_1)
    results_df_2 = extract_data_for_wilcoxon_and_t_test(cells[1], file_path_HE_2)
    
    cell_density_1_before = results_df_1["nbr_mm2_before"] + epsilon
    cell_density_2_before = results_df_2["nbr_mm2_before"] + epsilon
    
    ratio_before = cell_density_1_before/cell_density_2_before
    
    cell_density_1_after = results_df_1["nbr_mm2_after"] + epsilon
    cell_density_2_after = results_df_2["nbr_mm2_after"] + epsilon
    
    ratio_after = cell_density_1_after/cell_density_2_after
    
    # X positions for the bars
    x_positions = [1, 1.2]
    
    plt.style.use({'figure.facecolor': to_rgba('orchid', alpha=0.1), 'axes.facecolor': to_rgba('orchid', alpha=0.06)})
    plt.figure(figsize=[10,8])

    # # Bars
    # plt.bar(x_positions, max_vals, color=['darkorange', 'darkblue'], alpha=0.6, width=0.6)
    box = plt.boxplot(
                [ratio_before, ratio_after],
                positions=x_positions,
                patch_artist=True,    # fill the box with color
                medianprops=dict(color='black', linewidth=2),
                # whiskerprops=dict(color='gray', linewidth=1.5),
                capprops=dict(linewidth=0)  # removes the little line at whisker ends
            )
    
    # Box colors
    box_colors = ["darkorange", "maroon"]
    for patch, color in zip(box["boxes"], box_colors):
        patch.set_facecolor(color)
    # Overlay dot points
    # Add small random jitter so dots don’t overlap exactly
    jitter_before = np.random.normal(x_positions[0], 0.04, len(ratio_before))
    jitter_after = np.random.normal(x_positions[1], 0.04, len(ratio_after))
    x_pos_array_before = np.full(len(ratio_before), x_positions[0])
    x_pos_array_after = np.full(len(ratio_after), x_positions[1])
    
    plt.scatter(x_pos_array_before, ratio_before, alpha=0.6, color='black', s=20, zorder=3)
    plt.scatter(x_pos_array_after, ratio_after, alpha=0.6, color='black', s=20, zorder=3)
    
    save_path = Path(f"/media/jenny/Expansion/MM_HE_results/correlation/cell_ratios/all_biopsies/{cell}/")
    
    if not save_path.exists():
        save_path.mkdir(parents=True)
        print(f"Directory {save_path} was created")

    # Add labels
    # plt.text(f"P = {p:.1e}", fontsize=21, color="black", bbox=dict(facecolor="orchid",alpha=0.1, edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
    plt.xticks(x_positions, ["Before", "During"], fontsize=21)
    plt.yticks(fontsize=21)
    plt.ylabel(f"{cells[0]}-to-{cells[1]} ratio", fontsize=21, labelpad=12)
    plt.xlim(0.85, 1.35)
    plt.title(f"HoVer-Net", fontsize=21)
    plt.tight_layout()
    plt.savefig(save_path.joinpath(f"ratio_{cells[0]}_{cells[1]}.png"))
    plt.show()

if __name__ == "__main__":
    cell = "epithelial"
    file_path_HE = f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/all_statistics/all_statistics_results_{cell}_HE.csv"
    file_path_GT = f"/media/jenny/Expansion/MetoxyLacc_HE_20x_results/metoxylacc_tumor_percentage_heidi.csv"
    
    # file_path_HE = f"/media/jenny/Expansion1/MM_HE_results/all_statistics/excluding_discarded_statistics_results_{cell}_HE.csv"
    # file_path_CD8 = "/media/jenny/Expansion1/CellPopEstimatesIHC_QuPathAnalyses.csv"
    # file_path_genes = "/media/jenny/Expansion1/CellPopEstimatesGenes.csv"
    
    # t_test(file_path_HE = file_path_HE, cell=cell, column_list=["percent", "nbr_mm2", "tot_cells_per_mm2"], y_label_list=[f"% Epithelial cells", r"Epithelial cells per $\mu$m$^2$ (Log)", r"Cells per $\mu$m$^2$ (Log)"])
    t_test_gt_data(file_path_HE = file_path_GT, cell=cell, column_list=["percent"], y_label_list=[f"% Epithelial cells"])
    
    # wilcoxon_test(file_path_HE = file_path_HE, cell=cell, column_list=["percent", "nbr_mm2", "tot_cells_per_mm2"], y_label_list=[f"% Lymphocytes", r"Lymphocytes per $\mu$m$^2$ (Log)", r"Cells per $\mu$m$^2$ (Log)"])
    # pearson_test(patients=22, cell=cell, column_list=["T_cells_GeneEstimate", "ESTIMATE_ImmuneScore", "B_cells_GeneEstimate", "CD8_T_cells_GeneEstimate" , "Cytotoxic_cells_GeneEstimate", "NK_cells_GeneEstimate"], x_label_list =[f"% Lymphocytes",f"Lymphocytes per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # spearman_test(patients=22, cell=cell, column_list=["T_cells_GeneEstimate", "ESTIMATE_ImmuneScore", "B_cells_GeneEstimate", "CD8_T_cells_GeneEstimate" , "Cytotoxic_cells_GeneEstimate", "NK_cells_GeneEstimate"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # pearson_test(patients=22, cell=cell, column_list=["ESTIMATE_StromalScore"], x_label_list =[f"% Epithelial cells (IHC)",f"Epithelial cells per mm$^2$ (Augmented HoVer-Net)"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # spearman_test(patients=22, cell=cell, column_list=["ESTIMATE_StromalScore"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # pearson_test_gt_data(patients=40, cell=cell, column_list=["percentage_tumor_Heide"], x_label_list =[f"% Epithelial"], HE_column_list=[f"percent_nuclei_{cell}_cells"], file_path_HE=file_path_HE, file_path_GT=file_path_GT)
    
    # before_CD8(cell=cell, column_list=["aSMAPosCellsPercent_WholeTissue","NumberOfaSMAPosCellsPermm2_WholeTissue"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
    # after_CD8(cell=cell, column_list=["aSMAPosCellsPercent_WholeTissue","NumberOfaSMAPosCellsPermm2_WholeTissue"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
    # all_immune_cells_pearson_spearman(patients= 22, column_list=["ESTIMATE_ImmuneScore"], file_path_genes=file_path_genes)
    # all_immune_cells_pearson_spearman_plot_together(patients= 22, column_list=["ESTIMATE_ImmuneScore"], file_path_genes=file_path_genes)
    # tables()
    # before_and_after_CD8(cell=cell, column_list=["aSMAPosCellsPercent_WholeTissue","NumberOfaSMAPosCellsPermm2_WholeTissue"], y_label_list =[r'% $\alpha$SMA positive cells (IHC)', r'$\alpha$SMA cells per mm$^2 (IHC)$'], x_label_list =[f"% Stromal cells (Augmented HoVer-Net)",f"Stromal cells per mm$^2$ (Augmented HoVer-Net)"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
    # before_and_after_CD8(cell=cell, column_list=["NumberOfCD8PosCellsPermm2_WholeTissue", "NumberOfCD8PosCells_WholeTissue", "CD8PosCellsPercent_WholeTissue", "CD8PosPixelsPercent_WholeTissue"], y_label_list =['CD8 positive cells per mm$^2$', 'CD8 positive cells', '% CD8 positive cells', '% CD8 positive pixels', 'Cells per mm$^2$ (IHC)'], x_label_list =["Lymphocytes per mm$^2$", "Total number of lymphocytes", "% Lymphocytes", "% lymphocyte pixels"], HE_column_list=[f"nbr_{cell}_cells_per_mm2","nbr_nuclei_lymphocyte_cells", f"percent_nuclei_{cell}_cells", "percent_pixels_lymphocyte_cells"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
    # before_and_after_CD8(cell=cell, column_list=["NumberOfCellsPermm2_CD8section_WholeTissue"], y_label_list =['Cells per mm$^2$ (IHC)'], x_label_list =[ "Cells per mm$^2$ (Augmented HoVer-Net)"], HE_column_list=["tot_cells_per_mm2"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
   
    # before_and_after_CD8(cell=cell, column_list=["NumberOfCells_CD8section_WholeTissue"], y_label_list =['Total number of cells (IHC)'], x_label_list =["Total number of cells (Augmented HoVer-Net)"], HE_column_list=["tot_nuclei_wsi"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
    # cell_ratio(["neutrophil", "lymphocyte"])