import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


def read_table(path):
    return pd.read_csv(path, sep=None, engine="python")


def make_bounds(value, info):
    info = "" if pd.isna(info) else str(info).strip().lower()

    if info == "" or info == "exact":
        return float(value), float(value)
    elif info == "over":
        return float(value), 100.0
    elif info == "under":
        return 0.0, float(value)
    else:
        raise ValueError(f"Unknown additional_info value: {info}")


def interval_error(pred, lower, upper):
    if pred < lower:
        return lower - pred
    elif pred > upper:
        return pred - upper
    else:
        return 0.0


def effective_ground_truth_for_plot(pred, lower, upper):
    """
    For the 'interval is correct' plot:
    - exact labels stay exact
    - if prediction is inside the interval, use prediction as x-value
    - otherwise use the nearest interval boundary
    """
    if lower == upper:
        return lower
    elif pred < lower:
        return lower
    elif pred > upper:
        return upper
    else:
        return pred


def prepare_comparison(file1_path, file2_path, prediction_col):
    gt = read_table(file1_path)
    pred = read_table(file2_path)

    gt.columns = gt.columns.str.strip()
    pred.columns = pred.columns.str.strip()

    df = pd.merge(gt, pred, on="Biopsi", how="inner")

    bounds = df.apply(
        lambda row: make_bounds(row["percentage_tumor_Heide"], row["additional_info"]),
        axis=1
    )
    df[["gt_lower", "gt_upper"]] = pd.DataFrame(bounds.tolist(), index=df.index)

    df["prediction"] = df[prediction_col]

    df["interval_error"] = df.apply(
        lambda row: interval_error(row["prediction"], row["gt_lower"], row["gt_upper"]),
        axis=1
    )

    df["gt_for_interval_plot"] = df.apply(
        lambda row: effective_ground_truth_for_plot(
            row["prediction"], row["gt_lower"], row["gt_upper"]
        ),
        axis=1
    )

    exact_mask = df["gt_lower"] == df["gt_upper"]
    df["signed_error_exact_only"] = np.nan
    df.loc[exact_mask, "signed_error_exact_only"] = (
        df.loc[exact_mask, "prediction"] - df.loc[exact_mask, "gt_lower"]
    )

    return df


def summarize_results(df):
    exact_mask = df["gt_lower"] == df["gt_upper"]

    print("\nMerged data:")
    print(df[[
        "Biopsi",
        "percentage_tumor_Heide",
        "additional_info",
        "gt_lower",
        "gt_upper",
        "prediction",
        "interval_error"
    ]])

    print("\nInterval-based evaluation on all samples:")
    print(f"Mean interval error: {df['interval_error'].mean():.3f}")
    print(f"Median interval error: {df['interval_error'].median():.3f}")
    print(f"Max interval error: {df['interval_error'].max():.3f}")
    print(f"Fraction within allowed interval: {(df['interval_error'] == 0).mean():.3f}")

    if exact_mask.any():
        abs_err = np.abs(df.loc[exact_mask, "signed_error_exact_only"])
        rmse = np.sqrt(np.mean(df.loc[exact_mask, "signed_error_exact_only"] ** 2))
        bias = df.loc[exact_mask, "signed_error_exact_only"].mean()

        print("\nExact-ground-truth evaluation only:")
        print(f"Number of exact samples: {exact_mask.sum()}")
        print(f"MAE: {abs_err.mean():.3f}")
        print(f"RMSE: {rmse:.3f}")
        print(f"Bias: {bias:.3f}")


def plot_comparison_with_arrows(df, save_path_base, title="HoVerNet vs ground truth tumor percentage"):
    
    save_path = os.path.join(save_path_base, f"hovernet_gt_arrows.png")
    
    fig, ax = plt.subplots(figsize=(10, 13))

    for _, row in df.iterrows():
        pred = row["prediction"]
        value = row["percentage_tumor_Heide"]
        info = "" if pd.isna(row["additional_info"]) else str(row["additional_info"]).strip().lower()

        line, = ax.plot(value, pred, "o")
        color = line.get_color()

        if info == "over":
            ax.annotate(
                "",
                xy=(100, pred),
                xytext=(value, pred),
                arrowprops=dict(arrowstyle="->", lw=1.2, color=color)
            )
        elif info == "under":
            ax.annotate(
                "",
                xy=(0, pred),
                xytext=(value, pred),
                arrowprops=dict(arrowstyle="->", lw=1.2, color=color)
            )

        ax.text(value + 1, pred + 1, row["Biopsi"], fontsize=8, color=color)

    ax.plot([0, 100], [0, 100], "--")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Ground truth tumor percentage (interval-adjusted)")
    ax.set_ylabel("HoVerNet tumor percentage")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()


def plot_interval_as_correct(df, save_path_base, title="HoVerNet vs interval-adjusted ground truth"):
    
    save_path = os.path.join(save_path_base, f"hovernet_gt_threshold.png")
    
    fig, ax = plt.subplots(figsize=(10, 13))

    for _, row in df.iterrows():
        x = row["gt_for_interval_plot"]
        y = row["prediction"]

        line, = ax.plot(x, y, "o")
        color = line.get_color()
        ax.text(x + 1, y + 1, row["Biopsi"], fontsize=8, color=color)

    ax.plot([0, 100], [0, 100], "--")
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Ground truth tumor percentage (interval-adjusted)")
    ax.set_ylabel("HoVerNet tumor percentage")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()



if __name__ == "__main__":

    file2_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_results/all_statistics/all_statistics_results_epithelial_HE.csv"
    file1_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_results/metoxylacc_tumor_percentage_heidi.csv"
    
    prediction_col = "percent_nuclei_epithelial_cells"
    save_path_base = "/media/jenny/Expansion/MetoxyLacc_HE_20x_results/correlation/gt_comparison/"
    df = prepare_comparison(file1_path, file2_path, prediction_col)

    summarize_results(df)

    plot_comparison_with_arrows(
        df,
        save_path_base =save_path_base,
        title=f"Ground truth vs HoVerNet (with interval arrows)"
    )

    plot_interval_as_correct(
        df,
        save_path_base=save_path_base,
        title=f"Ground truth vs HoVerNet (interval treated as correct)"
    )

    df.to_csv("/media/jenny/Expansion/MetoxyLacc_HE_20x_results/correlation/gt_comparison/merged_comparison_results.csv", index=False)