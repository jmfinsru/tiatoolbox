import pandas as pd
import matplotlib.pyplot as plt


def analyze_amacr(input_path, output_path, output_path_amacr_positive):
    # Read data
    df = pd.read_csv(input_path)

    column_name = "Cell: AMACR OD mean"

    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' was not found in the file.")

    # Convert to numeric and remove missing/non-numeric values
    values = pd.to_numeric(
        df[column_name],
        errors="coerce"
    ).dropna()

    # -----------------------------------
    # Plot 1: All values
    # -----------------------------------
    plt.figure(figsize=(10, 6))

    plt.hist(
        values,
        bins=20,
        edgecolor="black"
    )

    plt.axvline(
        0,
        color="red",
        linestyle="--",
        label="Zero"
    )

    plt.xlabel("Cell: AMACR OD mean")
    plt.ylabel("Count")
    plt.title("Distribution of Cell: AMACR OD mean")
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()

    # -----------------------------------
    # Keep only values above 0.02
    # -----------------------------------
    amacr_above_002 = values[values > 0.02]

    if len(amacr_above_002) == 0:
        print("No values above 0.02 were found.")
        return

    # Calculate percentiles among values > 0.02
    amacr_95_percentile = amacr_above_002.quantile(0.95)
    amacr_99_percentile = amacr_above_002.quantile(0.99)

    print("-------------------------------")
    print(f"Number of values above 0.02: {len(amacr_above_002)}")
    print(f"95 percentile: {amacr_95_percentile}")
    print(f"99 percentile: {amacr_99_percentile}")
    print("-------------------------------")

    # -----------------------------------
    # Plot 2: Values above 0.02 only
    # -----------------------------------
    plt.figure(figsize=(10, 6))

    plt.hist(
        amacr_above_002,
        bins=20,
        edgecolor="black"
    )

    # Show percentile positions on plot
    plt.axvline(
        amacr_95_percentile,
        linestyle="--",
        label=f"95th percentile = {amacr_95_percentile:.4f}"
    )

    plt.axvline(
        amacr_99_percentile,
        linestyle="--",
        label=f"99th percentile = {amacr_99_percentile:.4f}"
    )

    plt.xlabel("Cell: AMACR OD mean")
    plt.ylabel("Count")
    plt.title("AMACR OD mean values above 0.02")
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_path_amacr_positive)
    # plt.show()

if __name__ == "__main__":
    
    biopsies = ["Func006","Func006_LN", "Func009", "Func015", "Func028", "Func029", "Func033", "Func036", "Func043", "Func044", "Func050", "Func050_LN", "Func052", "Func082", "Func083", "Func089", "Func097", "Func112", "Func116", "Func117"]
    # biopsies = ["Func117"]
    for biopsy in biopsies:
        print(biopsy)
        input_path = f"/media/jenny/Expansion/Prostata_QuPath/PIN_cells/{biopsy}/measurement_table/{biopsy}_PIN_trippel_Tumor_detections.csv"
        output_path =f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results_distribution_of_amacr/histogram_{biopsy}.png"
        output_path_amacr_positive = f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results_distribution_of_amacr/histogram_amacr_pos_{biopsy}.png"
        analyze_amacr(input_path, output_path, output_path_amacr_positive)