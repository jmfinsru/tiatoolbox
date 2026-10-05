import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from pathlib import Path
import tifffile


def load_spot_summary(summary_csv_path, column_name):
    """
    Load spot-level AMACR results.
    """
    df = pd.read_csv(summary_csv_path)

    required_columns = [
        "x_coreg",
        "y_coreg",
        f"{column_name}"
    ]

    missing = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )
    # Full-resolution TIFF dimensions

    return df


def get_spot_radius(radius_csv_path, biopsy):
    """
    Find the spot radius corresponding to biopsy
    in the nameList column.
    """
    radius_table = pd.read_csv(radius_csv_path)

    match = radius_table.loc[
        radius_table["nameList"] == biopsy,
        "spot_radius_coreg"
    ]

    if len(match) == 0:
        raise ValueError(
            f"{biopsy} was not found in the nameList column."
        )

    if len(match) > 1:
        raise ValueError(
            f"{biopsy} occurs more than once in the nameList column."
        )

    spot_radius = match.iloc[0]

    return float(spot_radius)


def load_downsampled_tif(
    tif_path,
    downsample_factor=10
):
    """
    Load a TIFF and downsample it for visualization.
    """
    print("Opening TIFF...")

    with tifffile.TiffFile(tif_path) as tif:
        page = tif.pages[0]

        print("Original TIFF shape:", page.shape)

        image = page.asarray()

    print("TIFF loaded.")

    if image.ndim == 2:

        image = image[
            ::downsample_factor,
            ::downsample_factor
        ]

    elif image.ndim == 3:

        # H x W x C
        if image.shape[-1] in [3, 4]:

            image = image[
                ::downsample_factor,
                ::downsample_factor,
                :
            ]

        # C x H x W
        elif image.shape[0] in [3, 4]:

            image = image[
                :,
                ::downsample_factor,
                ::downsample_factor
            ]

            image = np.moveaxis(
                image,
                0,
                -1
            )

        else:
            raise ValueError(
                f"Unsupported TIFF shape: {image.shape}"
            )

    else:
        raise ValueError(
            f"Unsupported TIFF dimensions: {image.shape}"
        )

    print("Downsampled image shape:", image.shape)

    return image


def create_overlay_image(
    tif_path,
    summary_csv_path,
    radius_csv_path,
    biopsy,
    output_image_path,
    column_name,
    downsample_factor=10,
    cmap_name="inferno",
    alpha=0.7
):
    """
    Create AMACR heatmap overlay on a TIFF image.

    The spot radius is automatically obtained from
    radius_csv_path using the biopsy name.
    """

    # --------------------------------------------------
    # Get spot radius
    # --------------------------------------------------

    spot_radius = get_spot_radius(
        radius_csv_path=radius_csv_path,
        biopsy=biopsy
    )

    print("Biopsy:", biopsy)
    print("Spot radius:", spot_radius)

    # --------------------------------------------------
    # Load AMACR spot results
    # --------------------------------------------------

    print("Loading spot data...")

    df = load_spot_summary(
        summary_csv_path, column_name
    )

    # --------------------------------------------------
    # Load TIFF
    # --------------------------------------------------

    image = load_downsampled_tif(
        tif_path,
        downsample_factor=downsample_factor
    )

    # --------------------------------------------------
    # Scale spot coordinates to match downsampled image
    # --------------------------------------------------

    df["display_x"] = (
        df["x_coreg"] /
        downsample_factor
    )

    df["display_y"] = (
        df["y_coreg"] /
        downsample_factor
    )

    display_radius = (
        spot_radius /
        downsample_factor
    )

    print(
        "Displayed spot radius:",
        display_radius
    )

    # --------------------------------------------------
    # Create figure
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(12, 12)
    )

    # --------------------------------------------------
    # Display TIFF
    # --------------------------------------------------

    if image.ndim == 2:

        ax.imshow(
            image,
            cmap="gray",
            origin="upper"
        )

    else:

        ax.imshow(
            image,
            origin="upper"
        )

    # --------------------------------------------------
    # Create spot circles
    # --------------------------------------------------

    patches = []

    for x, y in zip(
        df["display_x"],
        df["display_y"]
    ):

        circle = Circle(
            (x, y),
            radius=display_radius
        )

        patches.append(circle)

    # --------------------------------------------------
    # Heatmap values
    # --------------------------------------------------

    values = df[
        f"{column_name}"
    ].to_numpy()

    collection = PatchCollection(
        patches,
        cmap=cmap_name,
        alpha=alpha,
        linewidth=0
    )

    collection.set_array(values)

    collection.set_clim(
        0,
        1
    )

    ax.add_collection(
        collection
    )

    # --------------------------------------------------
    # Colorbar
    # --------------------------------------------------

    colorbar = fig.colorbar(
        collection,
        ax=ax,
        fraction=0.046,
        pad=0.04
    )

    colorbar.set_label(
        f"{column_name}"
    )

    # --------------------------------------------------
    # Image settings
    # --------------------------------------------------

    ax.set_xlim(
        0,
        image.shape[1]
    )

    ax.set_ylim(
        image.shape[0],
        0
    )

    ax.set_aspect(
        "equal"
    )

    ax.set_title(
        f"{biopsy}"
    )

    ax.axis(
        "off"
    )

    # --------------------------------------------------
    # Save overlay
    # --------------------------------------------------

    print("Saving overlay...")

    plt.savefig(
        output_image_path,
        dpi=200,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        "Overlay saved:",
        output_image_path
    )
    with tifffile.TiffFile(tif_path) as tif:
        original_shape = tif.pages[0].shape

        original_height = original_shape[0]
        original_width = original_shape[1]

        print("\n--- Coordinate diagnostic ---")

        print(
        "TIFF full-resolution size:",
        original_width,
        "x",
        original_height
        )

        print(
        "Spot X range:",
        df["x_coreg"].min(),
        "to",
        df["x_coreg"].max()
        )

        print(
        "Spot Y range:",
        df["y_coreg"].min(),
        "to",
        df["y_coreg"].max()
        )

        # Scale coordinates
        df["display_x"] = (
        df["x_coreg"] /
        downsample_factor
        )

        df["display_y"] = (
        df["y_coreg"] /
        downsample_factor
        )

        print(
        "Displayed X range:",
        df["display_x"].min(),
        "to",
        df["display_x"].max()
        )

        print(
        "Displayed Y range:",
        df["display_y"].min(),
        "to",
        df["display_y"].max()
        )

        print(
        "Downsampled image size:",
        image.shape[1],
        "x",
        image.shape[0]
        )

        print("-----------------------------\n")

    return df


if __name__ == "__main__":
    
    region_mode = "square" # circle or square, where square includes more tissue around the visium spots while the circle only look at the visium spots
    measurement_type = "cells" # cells or pixels
    I_ref = "0.5"
    column_name = "severity_adjusted_AMACR_mean"
    downsample_factor = 10
    biopsies = ["Func006","Func006_LN", "Func009", "Func015", "Func028", "Func029", "Func033", "Func036", "Func043", "Func044", "Func050", "Func050_LN", "Func052", "Func082", "Func083", "Func089", "Func097", "Func112", "Func116", "Func117"]
    # biopsies = ["Func117"]
    
    if measurement_type == "cells":
        for biopsy in biopsies:
            
            image_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/{biopsy}_PIN_trippel.tif"
            )
            
            radius_csv_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/spot_radius/spot_radius_coreg.csv"
            )

            summary_csv_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results/results_amacr_cells_in_spots/{biopsy}/{region_mode}/spot_AMACR_results_I_ref_{I_ref}.csv"
            )
            output_image_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results/results_amacr_cells_in_spots/{biopsy}/{region_mode}/{column_name}_overlay_heatmap_I_ref_{I_ref}.png"
            )

            # --------------------------------------------------
            # Run
            # --------------------------------------------------

            create_overlay_image(
                tif_path=image_path,
                summary_csv_path=summary_csv_path,
                radius_csv_path=radius_csv_path,
                biopsy=biopsy,
                output_image_path=output_image_path,
                downsample_factor=downsample_factor,
                column_name=column_name,
                cmap_name="jet",
                alpha=0.7
            )
    
    if measurement_type == "pixels":
        for biopsy in biopsies:
            
            image_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/{biopsy}_PIN_trippel.tif"
            )
            
            radius_csv_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/spot_radius/spot_radius_coreg.csv"
            )

            summary_csv_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results/results_amacr_pixel_in_spots/{biopsy}/{region_mode}/AMACR_pixel_measurements_I_ref_{I_ref}.csv"
            )
            output_image_path = Path(
                f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results/results_amacr_pixel_in_spots/{biopsy}/{region_mode}/{column_name}_overlay_heatmap_I_ref_{I_ref}.png"
            )

            # --------------------------------------------------
            # Run
            # --------------------------------------------------

            create_overlay_image(
                tif_path=image_path,
                summary_csv_path=summary_csv_path,
                radius_csv_path=radius_csv_path,
                biopsy=biopsy,
                output_image_path=output_image_path,
                downsample_factor=downsample_factor,
                column_name=column_name,
                cmap_name="jet",
                alpha=0.7
            )