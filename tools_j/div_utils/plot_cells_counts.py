from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def read_csv_auto_separator(path: Path) -> pd.DataFrame:
    """
    Reads either comma-separated or tab-separated CSV files.
    """
    return pd.read_csv(path, sep=None, engine="python")


def plot_cell_counts(INPUT_CSV, OUTPUT_FOLDER, OUTPUT_SUMMARY_CSV, OUTPUT_PLOT, CELL_TYPE_COLUMNS):
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    df = read_csv_auto_separator(INPUT_CSV)

    missing_cols = [col for col in CELL_TYPE_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in CSV: {missing_cols}")

    # Sum each cell-type column
    totals = df[CELL_TYPE_COLUMNS].sum()

    # Save totals as CSV
    totals_df = totals.reset_index()
    totals_df.columns = ["cell_type", "total_count"]
    totals_df.to_csv(OUTPUT_SUMMARY_CSV, index=False)

    print(f"Saved summary CSV to: {OUTPUT_SUMMARY_CSV}")

    # Remove zero-count categories from pie chart only
    pie_totals = totals[totals > 0]

    PINK_COLORS = [
    "#f7b6d2",  
    "#f48fb1",
    "#BC5B7B",
    "#d81b60",
    "#ad1457",  
    "#5b0027",
    ]

    pie_colors = [
        PINK_COLORS[CELL_TYPE_COLUMNS.index(cell_type)]
        for cell_type in pie_totals.index
    ]

    # Create plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Bar chart
    bars = axes[0].bar(
        totals.index,
        totals.values,
        color=PINK_COLORS,
        edgecolor="black",
        linewidth=0.5,
    )

    axes[0].set_title("Total cell counts per cell type")
    axes[0].set_ylabel("Cell counts")
    axes[0].tick_params(axis="x", rotation=30)

    # Increase y-axis limit so labels above bars are visible
    max_count = totals.max()
    if max_count > 0:
        axes[0].set_ylim(0, max_count * 1.20)
    else:
        axes[0].set_ylim(0, 1)

    # Add values above bars
    for bar, value in zip(bars, totals.values):
        axes[0].text(
            bar.get_x() + bar.get_width() / 2,
            value + max_count * 0.03 if max_count > 0 else 0.03,
            str(int(value)),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    # Pie chart
    if pie_totals.sum() > 0:
        axes[1].pie(
            pie_totals.values,
            labels=pie_totals.index,
            autopct="%1.0f %%",
            startangle=90,
            colors=pie_colors,
            textprops={"fontsize": 9},
        )
        axes[1].set_title("Cell type proportions")
    else:
        axes[1].text(
            0.5,
            0.5,
            "No non-zero cell counts",
            ha="center",
            va="center",
        )
        axes[1].set_title("Cell type proportions")

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved plot to: {OUTPUT_PLOT}")


if __name__ == "__main__":
    # Change these paths
    INPUT_CSV = Path("/media/jenny/Expansion1/jenny_funcprost/conic/results/Func116_ST_HE_20x_BF_01/counts/nuclei_counts_from_0_to_271.csv")
    OUTPUT_FOLDER = Path("/media/jenny/Expansion1/jenny_funcprost/conic/results/Func116_ST_HE_20x_BF_01/plots/")

    OUTPUT_SUMMARY_CSV = OUTPUT_FOLDER / "cell_type_totals.csv"
    OUTPUT_PLOT = OUTPUT_FOLDER / "cell_type_summary_plot.png"


    # CELL_TYPE_COLUMNS_PANNUKE = [
    #     "neoplastic",
    #     "inflammatory",
    #     "connective",
    #     "dead",
    #     "epithelial",
    # ]

    CELL_TYPE_COLUMNS_CONIC = [
        "neutrophil",
        "epithelial",
        "lymphocyte",
        "plasma",
        "eosinophil",
        "connective",
    ]

    plot_cell_counts(INPUT_CSV, OUTPUT_FOLDER, OUTPUT_SUMMARY_CSV, OUTPUT_PLOT, CELL_TYPE_COLUMNS_CONIC)