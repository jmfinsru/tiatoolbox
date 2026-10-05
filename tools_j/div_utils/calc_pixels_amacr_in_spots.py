import os
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


# ============================================================
# Load Visium spot information
# ============================================================

def load_spot_data(spots_path):
    """
    Load Visium spot positions and keep spots within tissue.
    """

    spots = pd.read_csv(spots_path)

    if "in_tissue" not in spots.columns:
        raise ValueError(
            "The spot file must contain an 'in_tissue' column."
        )

    spots = spots.loc[
        spots["in_tissue"] == 1
    ].copy()

    required = [
        "x_coreg",
        "y_coreg"
    ]

    missing = [
        col for col in required
        if col not in spots.columns
    ]

    if missing:
        raise ValueError(
            f"Missing spot coordinate columns: {missing}"
        )

    return spots


# ============================================================
# Load spot radius
# ============================================================

def get_spot_radius(
    radius_csv_path,
    biopsy
):
    """
    Get co-registered Visium spot radius.
    """

    radius_table = pd.read_csv(
        radius_csv_path
    )

    match = radius_table.loc[
        radius_table["nameList"] == biopsy,
        "spot_radius_coreg"
    ]

    if len(match) == 0:
        raise ValueError(
            f"{biopsy} was not found in nameList."
        )

    if len(match) > 1:
        raise ValueError(
            f"{biopsy} occurs more than once in nameList."
        )

    return float(
        match.iloc[0]
    )


# ============================================================
# Load AMACR TIFF
# ============================================================

def load_amacr_image(
    amacr_tif_path
):
    """
    Load the exported single-channel AMACR OD image.
    """

    print(
        f"Loading AMACR TIFF:\n{amacr_tif_path}"
    )

    image = tifffile.imread(
        amacr_tif_path
    )

    print(
        "Raw TIFF shape:",
        image.shape
    )

    # --------------------------------------------------------
    # Handle common single-channel layouts
    # --------------------------------------------------------

    if image.ndim == 2:

        amacr = image

    elif image.ndim == 3:

        # Possible:
        # 1 x H x W
        if image.shape[0] == 1:

            amacr = image[0]

        # H x W x 1
        elif image.shape[-1] == 1:

            amacr = image[..., 0]

        else:

            raise ValueError(
                f"Expected a single-channel AMACR image, "
                f"but got shape {image.shape}"
            )

    else:

        raise ValueError(
            f"Unsupported AMACR TIFF shape: {image.shape}"
        )

    amacr = amacr.astype(
        np.float32
    )

    print(
        "AMACR image shape:",
        amacr.shape
    )

    print(
        "AMACR minimum:",
        np.nanmin(amacr)
    )

    print(
        "AMACR maximum:",
        np.nanmax(amacr)
    )

    print(
        "AMACR mean:",
        np.nanmean(amacr)
    )

    return amacr


# ============================================================
# Pixel-level AMACR calculation
# ============================================================

def calculate_pixel_amacr_per_spot(
    spots,
    amacr_image,
    spot_radius,
    amacr_positive_threshold,
    amacr_i_ref,
    region_mode="square",
    rectangle_scale=1.3
):
    """
    Calculate pixel-level AMACR measurements for each Visium spot.

    region_mode:
        "square"
            MATLAB-style square:
            half_side = spot_radius * rectangle_scale

        "circle"
            Only pixels inside the actual Visium spot:
            (x - cx)^2 + (y - cy)^2 <= spot_radius^2

    Positive pixel:
        AMACR OD >= amacr_positive_threshold

    Severity-adjusted AMACR:
        positive_pixel_fraction *
        normalized_mean_positive_AMACR_OD
    """

    if amacr_i_ref <= amacr_positive_threshold:
        raise ValueError(
            "AMACR_I_REF must be greater than "
            "AMACR_POSITIVE_THRESHOLD."
        )

    region_mode = region_mode.lower()

    if region_mode not in ["square", "circle"]:
        raise ValueError(
            "region_mode must be either "
            "'square' or 'circle'."
        )

    image_height, image_width = (
        amacr_image.shape
    )

    spot_radius = float(spot_radius)

    # --------------------------------------------------------
    # Define extraction size
    # --------------------------------------------------------

    if region_mode == "square":

        half_side = round(
            spot_radius * float(rectangle_scale)
        )

        print("Region mode: SQUARE")

        print(
            f"Spot radius: {spot_radius:.3f} px"
        )

        print(
            f"Rectangle half-side: {half_side} px"
        )

        print(
            f"Rectangle size: "
            f"{2 * half_side + 1} x "
            f"{2 * half_side + 1} px"
        )

    else:

        # Bounding box around circle
        half_side = int(
            np.ceil(spot_radius)
        )

        print("Region mode: CIRCLE")

        print(
            f"Spot radius: {spot_radius:.3f} px"
        )

        print(
            f"Circle diameter: "
            f"{2 * spot_radius:.3f} px"
        )

    # --------------------------------------------------------
    # Prepare result arrays
    # --------------------------------------------------------

    results = spots.copy()

    n_spots = len(spots)

    pixel_count = np.zeros(
        n_spots,
        dtype=int
    )

    amacr_mean_all = np.full(
        n_spots,
        np.nan
    )

    amacr_sum_all = np.full(
        n_spots,
        np.nan
    )

    positive_pixel_count = np.zeros(
        n_spots,
        dtype=int
    )

    positive_pixel_fraction = np.full(
        n_spots,
        np.nan
    )

    positive_mean_od = np.full(
        n_spots,
        np.nan
    )

    positive_mean_od_normalized = np.full(
        n_spots,
        np.nan
    )

    positive_od_sum = np.full(
        n_spots,
        np.nan
    )

    severity_adjusted_amacr = np.full(
        n_spots,
        np.nan
    )

    # --------------------------------------------------------
    # Process each Visium spot
    # --------------------------------------------------------

    for j, (_, spot) in enumerate(
        spots.iterrows()
    ):

        cx = int(
            round(
                float(
                    spot["x_coreg"]
                )
            )
        )

        cy = int(
            round(
                float(
                    spot["y_coreg"]
                )
            )
        )

        # ----------------------------------------------------
        # Bounding box
        # ----------------------------------------------------

        x1 = max(
            cx - half_side,
            0
        )

        x2 = min(
            cx + half_side + 1,
            image_width
        )

        y1 = max(
            cy - half_side,
            0
        )

        y2 = min(
            cy + half_side + 1,
            image_height
        )

        spot_pixels = amacr_image[
            y1:y2,
            x1:x2
        ]

        # ----------------------------------------------------
        # Select either square or circle pixels
        # ----------------------------------------------------

        if region_mode == "square":

            # All pixels inside bounding rectangle
            values = spot_pixels.ravel()

        else:

            # ----------------------------------------------
            # Build coordinates corresponding to this crop
            # ----------------------------------------------

            yy, xx = np.ogrid[
                y1:y2,
                x1:x2
            ]

            # ----------------------------------------------
            # Circular mask centered on Visium spot
            # ----------------------------------------------

            circle_mask = (
                (xx - cx) ** 2
                +
                (yy - cy) ** 2
                <= spot_radius ** 2
            )

            # Keep only pixels inside circle
            values = spot_pixels[
                circle_mask
            ]

        # ----------------------------------------------------
        # Remove NaN / inf
        # ----------------------------------------------------

        values = values[
            np.isfinite(values)
        ]

        n_pixels = len(values)

        pixel_count[j] = n_pixels

        if n_pixels == 0:
            continue

        # ----------------------------------------------------
        # All-pixel measurements
        # ----------------------------------------------------

        amacr_mean_all[j] = np.mean(
            values
        )

        amacr_sum_all[j] = np.sum(
            values
        )

        # ----------------------------------------------------
        # AMACR-positive pixels
        # ----------------------------------------------------

        positive_mask = (
            values
            >= amacr_positive_threshold
        )

        positive_values = values[
            positive_mask
        ]

        n_positive = len(
            positive_values
        )

        positive_pixel_count[j] = (
            n_positive
        )

        fraction = (
            n_positive
            /
            n_pixels
        )

        positive_pixel_fraction[j] = (
            fraction
        )

        # ----------------------------------------------------
        # No AMACR-positive pixels
        # ----------------------------------------------------

        if n_positive == 0:

            positive_mean_od[j] = 0.0

            positive_mean_od_normalized[j] = 0.0

            positive_od_sum[j] = 0.0

            severity_adjusted_amacr[j] = 0.0

            continue

        # ----------------------------------------------------
        # Mean AMACR intensity of positive pixels
        # ----------------------------------------------------

        mean_positive = np.mean(
            positive_values
        )

        positive_mean_od[j] = (
            mean_positive
        )

        # ----------------------------------------------------
        # Normalize positive mean intensity
        #
        # This version matches the original MATLAB method:
        # lower bounded by 0, but NOT capped at 1.
        # ----------------------------------------------------

        mean_positive_norm = max(
            0.0,
            (
                mean_positive
                - amacr_positive_threshold
            )
            /
            (
                amacr_i_ref
                - amacr_positive_threshold
            )
        )

        positive_mean_od_normalized[j] = (
            mean_positive_norm
        )

        # ----------------------------------------------------
        # Sum of AMACR OD among positive pixels
        # ----------------------------------------------------

        positive_od_sum[j] = np.sum(
            positive_values
        )

        # ----------------------------------------------------
        # Severity-adjusted AMACR
        # ----------------------------------------------------

        severity_adjusted_amacr[j] = (
            fraction
            *
            mean_positive_norm
        )

    # --------------------------------------------------------
    # Add output columns
    # --------------------------------------------------------

    results[
        "measurement_region"
    ] = region_mode

    results[
        "spot_radius_px"
    ] = spot_radius

    if region_mode == "square":

        results[
            "rectangle_half_side_px"
        ] = half_side

    else:

        results[
            "rectangle_half_side_px"
        ] = np.nan

    results[
        "pixel_count"
    ] = pixel_count

    results[
        "AMACR_pixel_mean_all"
    ] = amacr_mean_all

    results[
        "AMACR_pixel_sum_all"
    ] = amacr_sum_all

    results[
        "AMACR_positive_pixel_count"
    ] = positive_pixel_count

    results[
        "AMACR_positive_pixel_fraction"
    ] = positive_pixel_fraction

    results[
        "AMACR_positive_pixel_mean_OD"
    ] = positive_mean_od

    results[
        "AMACR_positive_pixel_mean_OD_normalized"
    ] = positive_mean_od_normalized

    results[
        "AMACR_positive_pixel_OD_sum"
    ] = positive_od_sum

    results[
        "severity_adjusted_AMACR_pixel"
    ] = severity_adjusted_amacr

    return results
if __name__ == "__main__":

    # --------------------------------------------------------
    # AMACR settings
    # --------------------------------------------------------

    AMACR_POSITIVE_THRESHOLD = 0.02

    I_ref = 1.0

    REGION_MODE = "circle" # circle or square, where square includes more tissue around the visium spots while the circle only look at the visium spots

    RECTANGLE_SCALE = 1.3

    # --------------------------------------------------------
    # Shared radius file
    # --------------------------------------------------------

    RADIUS_CSV_PATH = (
        "/media/jenny/Expansion/"
        "Prostata_Vilde/Co-reg_images_10x/"
        "spot_radius/spot_radius_coreg.csv"
    )
    biopsies = ["Func006","Func006_LN", "Func009", "Func015", "Func028", "Func029", "Func033", "Func036", "Func043", "Func044", "Func050", "Func050_LN", "Func052", "Func082", "Func083", "Func089", "Func097", "Func112", "Func116", "Func117"]
    # biopsies = ["Func117"]
    for biopsy in biopsies:

        print()
        print("=" * 60)
        print(f"Processing {biopsy}")
        print("=" * 60)

        # ----------------------------------------------------
        # Input files
        # ----------------------------------------------------

        SPOTS_PATH = (
            "/media/jenny/Expansion/"
            "Prostata_Vilde/Co-reg_images_10x/"
            "Visium_spot_coordinates/"
            f"{biopsy}_Spot_Coordinates_CoReg.csv"
        )

        AMACR_TIF_PATH = (
            f"/media/jenny/Expansion/Prostata_QuPath/PIN_cells/{biopsy}/deconvolved_images/AMACR_deconvolved.ome.tif"
        )

        # ----------------------------------------------------
        # Output
        # ----------------------------------------------------

        OUTPUT_PATH = (
            "/media/jenny/Expansion/"
            "Prostata_Vilde/Co-reg_images_10x/"
            "results_amacr_pixel_in_spots/"
            f"{biopsy}/{REGION_MODE}/"
            f"AMACR_pixel_measurements_I_ref_{I_ref}.csv"
        )

        output_folder = os.path.dirname(
            OUTPUT_PATH
        )

        os.makedirs(
            output_folder,
            exist_ok=True
        )

        # ----------------------------------------------------
        # Load data
        # ----------------------------------------------------

        spots = load_spot_data(
            SPOTS_PATH
        )

        spot_radius = get_spot_radius(
            RADIUS_CSV_PATH,
            biopsy
        )

        amacr_image = load_amacr_image(
            AMACR_TIF_PATH
        )

        # ----------------------------------------------------
        # Coordinate diagnostic
        # ----------------------------------------------------

        print()
        print("--- Coordinate diagnostic ---")

        print(
            "AMACR image width:",
            amacr_image.shape[1]
        )

        print(
            "AMACR image height:",
            amacr_image.shape[0]
        )

        print(
            "Spot X:",
            spots["x_coreg"].min(),
            "to",
            spots["x_coreg"].max()
        )

        print(
            "Spot Y:",
            spots["y_coreg"].min(),
            "to",
            spots["y_coreg"].max()
        )

        print("-----------------------------")
        print()

        # ----------------------------------------------------
        # Calculate pixel-level AMACR
        # ----------------------------------------------------

        results = calculate_pixel_amacr_per_spot(
            spots=spots,
            amacr_image=amacr_image,
            spot_radius=spot_radius,
            amacr_positive_threshold=AMACR_POSITIVE_THRESHOLD,
            amacr_i_ref=I_ref,
            region_mode=REGION_MODE,
            rectangle_scale=RECTANGLE_SCALE
        )

        # ----------------------------------------------------
        # Save
        # ----------------------------------------------------

        results.to_csv(
            OUTPUT_PATH,
            index=False
        )

        print(
            f"Saved:\n{OUTPUT_PATH}"
        )