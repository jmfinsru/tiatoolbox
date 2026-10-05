import numpy as np
import pandas as pd
import os
from scipy.spatial import cKDTree


# ============================================================
# Data loading
# ============================================================

def load_spot_data(spots_path):
    """
    Load Visium spot positions and retain only spots within tissue.
    """
    spots = pd.read_csv(spots_path)

    if "in_tissue" not in spots.columns:
        raise ValueError(
            "The spot file must contain an 'in_tissue' column."
        )

    spots = spots.loc[
        spots["in_tissue"] == 1
    ].copy()

    return spots


def get_spot_radius(radius_csv_path, biopsy):
    """
    Find the spot radius corresponding to biopsy in the nameList column.
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

    spot_radius = float(match.iloc[0])

    return spot_radius


def load_amacr_cell_data(cell_csv_path):
    """
    Load the cell-level AMACR data.

    Required columns:
        Nucleus X px
        Nucleus Y px
        Cell: AMACR OD mean
    """
    required_columns = [
        "Nucleus X px",
        "Nucleus Y px",
        "Cell: AMACR OD mean",
    ]

    cells = pd.read_csv(
        cell_csv_path,
        usecols=required_columns
    )

    # Force required columns to numeric
    for column in required_columns:
        cells[column] = pd.to_numeric(
            cells[column],
            errors="coerce"
        )

    n_before = len(cells)

    # A cell cannot be used if coordinates or AMACR OD are missing
    cells = cells.dropna(
        subset=required_columns
    ).copy()

    n_removed = n_before - len(cells)

    if n_removed > 0:
        print(
            f"Removed {n_removed} cells with missing/non-numeric "
            f"coordinates or AMACR OD."
        )

    return cells


# ============================================================
# Spot coordinates
# ============================================================
def get_spot_coordinate_columns(spots):
    """
    Identify the x/y coordinate columns.

    Prefer co-registered coordinates when available.
    """

    possible_columns = [
        ("x_coreg", "y_coreg"),
        ("x_fullres", "y_fullres"),
        ("pxl_col_in_fullres", "pxl_row_in_fullres"),
        ("imagecol", "imagerow"),
    ]

    for x_col, y_col in possible_columns:
        if x_col in spots.columns and y_col in spots.columns:
            print(f"Using spot coordinates: {x_col}, {y_col}")
            return x_col, y_col

    raise ValueError(
        "Could not identify the spot coordinate columns.\n"
        "Expected one of:\n"
        "  x_coreg / y_coreg\n"
        "  x_fullres / y_fullres\n"
        "  pxl_col_in_fullres / pxl_row_in_fullres\n"
        "  imagecol / imagerow"
    )

def calculate_spot_amacr(
    spots,
    cells,
    spot_radius,
    amacr_positive_threshold,
    amacr_i_ref,
    region_mode="circle",
    radius_scale=1.0,
    rectangle_scale=1.3
):
    """
    Calculate AMACR measurements for each Visium spot.

    Parameters
    ----------
    region_mode : str
        "circle"
            Cells are included if the nucleus is within:

                spot_radius * radius_scale

            of the spot center.

        "square"
            Cells are included inside a MATLAB-style square:

                half_side = round(
                    spot_radius * rectangle_scale
                )

            With rectangle_scale=1.3, this reproduces:

                spot_radius + 0.3 * spot_radius

    Normalization
    -------------
    Positive AMACR values are normalized as:

        normalized_OD =
            max(
                0,
                (OD - threshold)
                /
                (I_ref - threshold)
            )

    There is NO upper clipping.

    Therefore values can be > 1 when OD > I_ref.
    """

    # --------------------------------------------------------
    # Check settings
    # --------------------------------------------------------

    if amacr_i_ref <= amacr_positive_threshold:
        raise ValueError(
            "amacr_i_ref must be greater than "
            "amacr_positive_threshold."
        )

    region_mode = region_mode.lower()

    if region_mode not in ["circle", "square"]:
        raise ValueError(
            "region_mode must be either "
            "'circle' or 'square'."
        )

    # --------------------------------------------------------
    # Spot coordinate columns
    # --------------------------------------------------------

    x_col, y_col = get_spot_coordinate_columns(
        spots
    )

    spot_radius = float(
        spot_radius
    )

    # --------------------------------------------------------
    # Cell coordinates
    # --------------------------------------------------------

    cell_coordinates = cells[
        [
            "Nucleus X px",
            "Nucleus Y px"
        ]
    ].to_numpy(
        dtype=float
    )

    cell_x = cell_coordinates[:, 0]
    cell_y = cell_coordinates[:, 1]

    cell_tree = cKDTree(
        cell_coordinates
    )

    amacr_values = cells[
        "Cell: AMACR OD mean"
    ].to_numpy(
        dtype=float
    )

    # --------------------------------------------------------
    # Define measurement geometry
    # --------------------------------------------------------

    if region_mode == "circle":

        measurement_radius = (
            spot_radius
            * float(radius_scale)
        )

        print("Measurement region: CIRCLE")
        print(
            f"Spot radius: "
            f"{spot_radius:.3f} px"
        )
        print(
            f"Radius scale: "
            f"{radius_scale}"
        )
        print(
            f"Measurement radius: "
            f"{measurement_radius:.3f} px"
        )

        half_side = np.nan

    else:

        half_side = round(
            spot_radius
            * float(rectangle_scale)
        )

        print("Measurement region: SQUARE")
        print(
            f"Spot radius: "
            f"{spot_radius:.3f} px"
        )
        print(
            f"Rectangle scale: "
            f"{rectangle_scale}"
        )
        print(
            f"Rectangle half-side: "
            f"{half_side} px"
        )
        print(
            f"Rectangle size: "
            f"{2 * half_side + 1} x "
            f"{2 * half_side + 1} px"
        )

    # --------------------------------------------------------
    # Output arrays
    # --------------------------------------------------------

    n_spots = len(
        spots
    )

    total_cells = np.zeros(
        n_spots,
        dtype=int
    )

    positive_cells = np.zeros(
        n_spots,
        dtype=int
    )

    positive_fraction = np.full(
        n_spots,
        np.nan
    )

    mean_positive_od_raw = np.full(
        n_spots,
        np.nan
    )

    mean_positive_od_norm = np.full(
        n_spots,
        np.nan
    )

    sum_positive_od_raw = np.full(
        n_spots,
        np.nan
    )

    sum_positive_od_norm = np.full(
        n_spots,
        np.nan
    )

    severity_adjusted_mean = np.full(
        n_spots,
        np.nan
    )

    severity_adjusted_sum = np.full(
        n_spots,
        np.nan
    )

    # --------------------------------------------------------
    # Process spots
    # --------------------------------------------------------

    for j, (_, spot) in enumerate(
        spots.iterrows()
    ):

        cx = float(
            spot[x_col]
        )

        cy = float(
            spot[y_col]
        )

        # ====================================================
        # CIRCLE
        # ====================================================

        if region_mode == "circle":

            cell_indices = (
                cell_tree.query_ball_point(
                    [cx, cy],
                    r=measurement_radius
                )
            )

            cell_indices = np.asarray(
                cell_indices,
                dtype=int
            )

        # ====================================================
        # SQUARE
        # ====================================================

        else:

            # Search circle that contains the whole square
            candidate_radius = (
                np.sqrt(2)
                * half_side
            )

            candidate_indices = (
                cell_tree.query_ball_point(
                    [cx, cy],
                    r=candidate_radius
                )
            )

            candidate_indices = np.asarray(
                candidate_indices,
                dtype=int
            )

            if len(candidate_indices) > 0:

                candidate_x = cell_x[
                    candidate_indices
                ]

                candidate_y = cell_y[
                    candidate_indices
                ]

                inside_square = (
                    (
                        np.abs(
                            candidate_x - cx
                        )
                        <= half_side
                    )
                    &
                    (
                        np.abs(
                            candidate_y - cy
                        )
                        <= half_side
                    )
                )

                cell_indices = (
                    candidate_indices[
                        inside_square
                    ]
                )

            else:

                cell_indices = np.array(
                    [],
                    dtype=int
                )

        # ----------------------------------------------------
        # Total cells
        # ----------------------------------------------------

        n_total = len(
            cell_indices
        )

        total_cells[j] = (
            n_total
        )

        if n_total == 0:

            positive_cells[j] = 0
            positive_fraction[j] = np.nan

            mean_positive_od_raw[j] = np.nan
            mean_positive_od_norm[j] = np.nan

            sum_positive_od_raw[j] = 0.0
            sum_positive_od_norm[j] = 0.0

            severity_adjusted_mean[j] = np.nan
            severity_adjusted_sum[j] = np.nan

            continue

        # ----------------------------------------------------
        # AMACR values inside region
        # ----------------------------------------------------

        spot_amacr = amacr_values[
            cell_indices
        ]

        # ----------------------------------------------------
        # AMACR-positive cells
        # ----------------------------------------------------

        positive_mask = (
            spot_amacr
            >= amacr_positive_threshold
        )

        positive_raw_values = (
            spot_amacr[
                positive_mask
            ]
        )

        n_positive = len(
            positive_raw_values
        )

        positive_cells[j] = (
            n_positive
        )

        # ----------------------------------------------------
        # Positive fraction
        # ----------------------------------------------------

        fraction = (
            n_positive
            /
            n_total
        )

        positive_fraction[j] = (
            fraction
        )

        # ----------------------------------------------------
        # No AMACR-positive cells
        # ----------------------------------------------------

        if n_positive == 0:

            mean_positive_od_raw[j] = 0.0
            mean_positive_od_norm[j] = 0.0

            sum_positive_od_raw[j] = 0.0
            sum_positive_od_norm[j] = 0.0

            severity_adjusted_mean[j] = 0.0
            severity_adjusted_sum[j] = 0.0

            continue

        # ----------------------------------------------------
        # Raw AMACR measurements
        # ----------------------------------------------------

        mean_raw = np.mean(
            positive_raw_values
        )

        sum_raw = np.sum(
            positive_raw_values
        )

        mean_positive_od_raw[j] = (
            mean_raw
        )

        sum_positive_od_raw[j] = (
            sum_raw
        )

        # ----------------------------------------------------
        # Normalize EACH positive AMACR value
        #
        # threshold -> 0
        # I_ref     -> 1
        #
        # IMPORTANT:
        # No upper clipping.
        #
        # Values > I_ref therefore give normalized values > 1.
        # ----------------------------------------------------

        positive_normalized_values = (
            positive_raw_values
            - amacr_positive_threshold
        ) / (
            amacr_i_ref
            - amacr_positive_threshold
        )

        # Only lower-bound at zero
        positive_normalized_values = (
            np.maximum(
                positive_normalized_values,
                0
            )
        )

        # ----------------------------------------------------
        # Mean normalized AMACR OD
        # ----------------------------------------------------

        mean_norm = np.mean(
            positive_normalized_values
        )

        mean_positive_od_norm[j] = (
            mean_norm
        )

        # ----------------------------------------------------
        # Sum normalized AMACR OD
        # ----------------------------------------------------

        sum_norm = np.sum(
            positive_normalized_values
        )

        sum_positive_od_norm[j] = (
            sum_norm
        )

        # ----------------------------------------------------
        # Severity-adjusted mean AMACR
        #
        # Closest analogue to original SAHF:
        #
        # positive fraction
        # ×
        # normalized mean positive intensity
        # ----------------------------------------------------

        severity_adjusted_mean[j] = (
            fraction
            *
            mean_norm
        )

        # ----------------------------------------------------
        # Severity-adjusted sum AMACR
        # ----------------------------------------------------

        severity_adjusted_sum[j] = (
            fraction
            *
            sum_norm
        )

    # --------------------------------------------------------
    # Build output dataframe
    # --------------------------------------------------------

    results = spots.copy()

    results[
        "measurement_region"
    ] = region_mode

    results[
        "spot_radius_px"
    ] = spot_radius

    if region_mode == "circle":

        results[
            "measurement_radius_px"
        ] = measurement_radius

        results[
            "rectangle_half_side_px"
        ] = np.nan

    else:

        results[
            "measurement_radius_px"
        ] = np.nan

        results[
            "rectangle_half_side_px"
        ] = half_side

    results[
        "total_cells"
    ] = total_cells

    results[
        "AMACR_positive_cells"
    ] = positive_cells

    results[
        "AMACR_positive_fraction"
    ] = positive_fraction

    results[
        "AMACR_positive_mean_OD_raw"
    ] = mean_positive_od_raw

    results[
        "AMACR_positive_mean_OD_normalized"
    ] = mean_positive_od_norm

    results[
        "AMACR_positive_OD_sum_raw"
    ] = sum_positive_od_raw

    results[
        "AMACR_positive_OD_sum_normalized"
    ] = sum_positive_od_norm

    results[
        "severity_adjusted_AMACR_mean"
    ] = severity_adjusted_mean

    results[
        "severity_adjusted_AMACR_sum"
    ] = severity_adjusted_sum

    return results

if __name__ == "__main__":
    
    # Cell is considered AMACR positive at or above this OD
    AMACR_POSITIVE_THRESHOLD = 0.02 

    # Reference OD corresponding to normalized intensity = 1
    AMACR_I_REF = 1.0
    
    REGION_MODE = "square" # circle or square, where square includes more tissue around the visium spots while the circle only look at the visium spots

    RADIUS_SCALE = 1.0
    RECTANGLE_SCALE = 1.3

    biopsies = ["Func006","Func006_LN", "Func009", "Func015", "Func028", "Func029", "Func033", "Func036", "Func043", "Func044", "Func050", "Func050_LN", "Func052", "Func082", "Func083", "Func089", "Func097", "Func112", "Func116", "Func117"]
    # biopsies = ["Func117"]
    for biopsy in biopsies:

        SPOTS_PATH = f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/Visium_spot_coordinates/{biopsy}_Spot_Coordinates_CoReg.csv"
        RADIUS_CSV_PATH = "/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/spot_radius/spot_radius_coreg.csv"
        AMACR_CELL_CSV_PATH = f"/media/jenny/Expansion/Prostata_QuPath/PIN_cells/{biopsy}/measurement_table/{biopsy}_PIN_trippel_Tumor_detections.csv"
        OUTPUT_PATH = f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results/results_amacr_cells_in_spots/{biopsy}/{REGION_MODE}/spot_AMACR_results_I_ref_{AMACR_I_REF}.csv"


        spots = load_spot_data(
            SPOTS_PATH
        )

        cells = load_amacr_cell_data(
            AMACR_CELL_CSV_PATH
        )

        spot_radius = get_spot_radius(
            RADIUS_CSV_PATH,
            biopsy
        )

        results = calculate_spot_amacr(
            spots=spots,
            cells=cells,
            spot_radius=spot_radius,
            amacr_positive_threshold=AMACR_POSITIVE_THRESHOLD,
            amacr_i_ref=AMACR_I_REF,
            region_mode=REGION_MODE,
            radius_scale=RADIUS_SCALE,
            rectangle_scale=RECTANGLE_SCALE
        )

        # Create output folder if it does not exist
        output_folder = os.path.dirname(OUTPUT_PATH)

        if output_folder:
            os.makedirs(output_folder, exist_ok=True)

        # Save results
        results.to_csv(
            OUTPUT_PATH,
            index=False
        )

        print(f"Saved results to: {OUTPUT_PATH}")