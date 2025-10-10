import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd

from scipy import stats
from scipy.stats import wilcoxon, ttest_rel
from pathlib import Path

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
            percent = df_HE.loc[df_HE["Biopsi"] == biopsi, f"percent_nuclei_{cell}_cells"].values
            nbr_type_cells_per_mm2 = df_HE.loc[df_HE["Biopsi"] == biopsi, f"nbr_{cell}_cells_per_mm2"].values
            tot_cells_per_mm2 = df_HE.loc[df_HE["Biopsi"] == biopsi, "tot_cells_per_mm2"].values
            CD8_dict.update({
                    "Biopsi": biopsi,
                    "Group": group,
                    f"percent_nuclei_{cell}_cells": float(percent),
                    f"nbr_{cell}_cells_per_mm2" : float(nbr_type_cells_per_mm2), 
                    "tot_cells_per_mm2" : float(tot_cells_per_mm2),
                })

            results.append(CD8_dict.copy())

    # Create a dataframe with the results and reorder the columns
    results_df = pd.DataFrame(results)
    new_order_part_1 = ["Biopsi", "Group", f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2", "tot_cells_per_mm2" ]
    new_order = new_order_part_1 + column_list
    results_df = results_df[new_order]
    
    return results_df


def before_CD8_lymph(cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_CD8: str | Path):
    
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

def after_CD8_lymph( cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_CD8: str | Path ):
    
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


def extract_gene_data_spearman_and_pearson(column_list: list, file_path_HE: str | Path, file_path_genes: str | Path):
    
    # Load data files as a Dataframe
    df_HE = pd.read_csv(file_path_HE)
    df_genes = pd.read_csv(file_path_genes)

    # Regex patterns for the Biopsies
    pattern_before = re.compile(r"\.[A-Z]$")
    pattern_after = re.compile(r"\.\d+[A-Z]$")

    results = []
    # Loop through all biopsies in the genes Dataframe
    for gathered_biopsies in df_genes["Biopsi"]:
        genes_dict = {}
        for column in column_list: # Loop through columns of interest to extract the data
                    col = df_genes.loc[df_genes["Biopsi"]==gathered_biopsies, column].values
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
    

def pearson_test(cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_genes: str | Path):           
    
    results_df = extract_gene_data_spearman_and_pearson(column_list, file_path_HE, file_path_genes)
    # Loop through all column names in column_list
    for i in range(len(column_list)):
        for j in range(len(HE_column_list)):
            column = column_list[i]
            column_HE = HE_column_list[j]
            
            # Calculate pearson correlation for before and after treatment
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b,p_b = stats.pearsonr(before_treatment[column_HE], before_treatment[column])
            
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a,p_a = stats.pearsonr(after_treatment[column_HE], after_treatment[column])

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/pearson/genes_{cell}/{column_HE}")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")

            plt.figure(figsize=(12, 10))
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "rebeccapurple", label =f"Pearson  P={p_b:.3f}, r={r_b:.2f}")
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
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "rebeccapurple", label =f"Pearson  P={p_a:.3f}, r={r_a:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_pearson_{column}_after.png"))
            plt.show()

def spearman_test(cell: str, x_label_list: list, column_list: list, HE_column_list: list, file_path_HE: str | Path, file_path_genes: str | Path):           
    
    results_df = extract_gene_data_spearman_and_pearson(column_list, file_path_HE, file_path_genes)
    # Loop through all column names in column_list
    for i in range(len(column_list)):
        for j in range(len(HE_column_list)):
            column = column_list[i]
            column_HE = HE_column_list[j]
            
            # Calculate spearman correlation for before and after treatment
            before_treatment = results_df[results_df["Group"] == "before"]
            r_b,p_b = stats.spearmanr(before_treatment[column_HE], before_treatment[column])
            
            after_treatment = results_df[results_df["Group"] == "after"]
            r_a,p_a = stats.spearmanr(after_treatment[column_HE], after_treatment[column])

            # Path to save plots
            save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/spearman/genes_{cell}/{column_HE}")
            if not save_path.exists():
                save_path.mkdir(parents=True)
                print(f"Directory {save_path} was created")

            plt.figure(figsize=(12, 10))
            plt.scatter(before_treatment[column_HE], before_treatment[column], color = "rebeccapurple", label =f"Spearman  P={p_b:.3f}, r={r_b:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(before_treatment["Biopsi"],before_treatment[column_HE], before_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_spearman_{column}_before.png"))
            plt.show()

            plt.figure(figsize=(12, 10))
            plt.scatter(after_treatment[column_HE], after_treatment[column], color = "rebeccapurple", label =f"Spearman  P={p_a:.3f}, r={r_a:.2f}")
            plt.xlabel(x_label_list[j])
            plt.ylabel(column)
            # plt.text(7.5, 1750, f"P = {p1:.2f}, r = {r1:.2f}", fontsize=10, color="black", bbox=dict(facecolor="thistle", edgecolor="black", boxstyle="round,pad=0.5"), alpha = 0.9)
            plt.legend()

            # Add labels with Biopsi IDs
            for k, x, y in zip(after_treatment["Biopsi"],after_treatment[column_HE], after_treatment[column]):
                plt.text(x, y, k, fontsize=8)
            plt.savefig(save_path.joinpath(f"corr_spearman_{column}_after.png"))
            plt.show()


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
    print(results_df)
    
    # Loop through all column names in column_list
    for i in range(len(column_list)):
        
        column = column_list[i]
        y_label = y_label_list[i]
        before_treatment = np.concatenate(results_df[column + "_before"].values)
        after_treatment = np.concatenate(results_df[column + "_after"].values)
        
        # print("----------------------------------------------------------------------------")
        # print(f"before_treatment: {before_treatment}")
        # print(f"after_treatment: {after_treatment}")
        # print("----------------------------------------------------------------------------")
    
        # Calculate the Wilcoxon signed-rank test
        try:
            stat, p = wilcoxon(before_treatment, after_treatment)
        except ValueError:
            stat, p = None, None

        print(f"Wilcoxon signed-rank test on {y_label} gives: stat={stat},  p={p}")
        
        # X positions for the bars
        x_positions = [0, 1]

        # Heights of bars = max values
        max_vals = [before_treatment.max(), after_treatment.max()]
        
        plt.figure(figsize=[12,9])
        # Bars
        plt.bar(x_positions, max_vals, color=['skyblue', 'darkblue'], alpha=0.3, width=0.6)

        # Plot individual dots
        for i in range(len(before_treatment)):
            plt.plot(x_positions, [before_treatment[i], after_treatment[i]], marker='o', color='black', alpha=0.7)
            plt.text(-0.12, before_treatment[i], results_df["before"][i], fontsize=8)
            plt.text(1.02, after_treatment[i], results_df["after"][i], fontsize=8)
        # Path to save results
        save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/wilcoxon/{cell}")
        
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")

        # Add labels
        plt.xticks(x_positions, ["Before", "After"])
        plt.ylabel(y_label)
        plt.title(f"Wilcoxon test, p = {p:.6f}")
        plt.savefig(save_path.joinpath(f"wilcoxon_{column}.png"))
        plt.show()

def t_test(column_list: list, y_label_list: list, cell: str, file_path_HE: str | Path):
    
    # Load dataframe
    results_df = extract_data_for_wilcoxon_and_t_test(cell, file_path_HE)

    # Loop through all column names in column_list
    for i in range(len(column_list)):
        
        column = column_list[i]
        y_label = y_label_list[i]
        before_treatment = np.concatenate(results_df[column + "_before"].values)
        after_treatment = np.concatenate(results_df[column + "_after"].values)
       
        # print("----------------------------------------------------------------------------")
        # print(f"before_treatment: {before_treatment}")
        # print(f"after_treatment: {after_treatment}")
        # print("----------------------------------------------------------------------------")
    
        # Perform t-test
        try:
            stat, p = ttest_rel(before_treatment, after_treatment)
        except ValueError:
            stat, p = None, None

        print(f"t-test on {y_label} gives: stat={stat},  p={p}")
        
        # X positions for the bars
        x_positions = [0, 1]

        # Heights of bars = max values
        max_vals = [before_treatment.max(), after_treatment.max()]
        
        plt.figure(figsize=[12,9])
        plt.bar(x_positions, max_vals, color=['skyblue', 'darkblue'], alpha=0.3, width=0.6)

        # Plot individual dots
        for i in range(len(before_treatment)):
            plt.plot(x_positions, [before_treatment[i], after_treatment[i]], marker='o', color='black', alpha=0.7)
            plt.text(-0.12, before_treatment[i], results_df["before"][i], fontsize=8)
            plt.text(1.02, after_treatment[i], results_df["after"][i], fontsize=8)
        
        save_path = Path(f"/media/jenny/Expansion/MM_HE_results/Correlation/t_test/{cell}")
        if not save_path.exists():
            save_path.mkdir(parents=True)
            print(f"Directory {save_path} was created")

        # Add labels
        plt.xticks(x_positions, ["Before", "After"])
        plt.ylabel(y_label)
        plt.title(f"t-test, p = {p:.6f}")
        plt.savefig(save_path.joinpath(f"t_test_{column}.png"))
        plt.show()
        

if __name__ == "__main__":
    cell = "connective"
    file_path_HE = f"/media/jenny/Expansion/MM_HE_results/all_statistics/all_statistics_{cell}_HE.csv"
    file_path_CD8 = "/media/jenny/Expansion/CellPopEstimatesIHC_QuPathAnalyses.csv"
    file_path_genes = "/media/jenny/Expansion/CellPopEstimatesGenes.csv"
    # t_test(file_path_HE = file_path_HE, cell=cell, column_list=["percent", "nbr_mm2", "tot_cells_per_mm2"], y_label_list=[f"% {cell}", f"{cell} per mm$^2$", "cells per mm$^2$"])
    # wilcoxon_test(file_path_HE = file_path_HE, cell=cell, column_list=["percent", "nbr_mm2", "tot_cells_per_mm2"], y_label_list=[f"% {cell}", f"{cell} per mm$^2$", "cells per mm$^2$"])
    # pearson_test(cell=cell, column_list=["T_cells_GeneEstimate", "ESTIMATE_ImmuneScore", "B_cells_GeneEstimate", "CD8_T_cells_GeneEstimate" , "Cytotoxic_cells_GeneEstimate", "NK_cells_GeneEstimate"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # spearman_test(cell=cell, column_list=["T_cells_GeneEstimate", "ESTIMATE_ImmuneScore", "B_cells_GeneEstimate", "CD8_T_cells_GeneEstimate" , "Cytotoxic_cells_GeneEstimate", "NK_cells_GeneEstimate"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # pearson_test(cell=cell, column_list=["ESTIMATE_StromalScore"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # spearman_test(cell=cell, column_list=["ESTIMATE_StromalScore"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_genes=file_path_genes)
    # before_CD8_lymph(cell=cell, column_list=["aSMAPosCellsPercent_WholeTissue","NumberOfaSMAPosCellsPermm2_WholeTissue"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
    after_CD8_lymph(cell=cell, column_list=["aSMAPosCellsPercent_WholeTissue","NumberOfaSMAPosCellsPermm2_WholeTissue"], x_label_list =[f"% {cell} cells",f"number of {cell} cells per mm$^2$"], HE_column_list=[f"percent_nuclei_{cell}_cells", f"nbr_{cell}_cells_per_mm2"], file_path_HE=file_path_HE, file_path_CD8=file_path_CD8)
