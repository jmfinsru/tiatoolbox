import os
import re
from pathlib import Path

import cv2
import numpy as np
from natsort import natsorted


def extract_patch_number(filename):
    """
    Extract patch number from filenames such as:

        epithelial_patch_7.png
        neutrophil_patch_7.png
        overlay_patch_7.png

    Returns:
        7
    """

    match = re.search(r"patch_(\d+)", filename)

    if match is None:
        return None

    return int(match.group(1))


def remove_boundaries(
    filled_path,
    boundary_path,
    output_path,
    boundary_threshold=0,
):
    """
    Remove all colored/non-black boundary pixels from a binary filled
    nuclei mask.

    Filled mask:
        0   = background
        255 = nucleus

    Boundary image:
        black = no boundary
        anything non-black = boundary

    The output remains binary:
        0   = background / removed boundary
        255 = nucleus
    """

    # ---------------------------------------------------------
    # Read filled binary nuclei mask
    # ---------------------------------------------------------

    filled = cv2.imread(
        str(filled_path),
        cv2.IMREAD_GRAYSCALE
    )

    if filled is None:
        raise ValueError(
            f"Could not read filled nuclei image:\n{filled_path}"
        )

    # ---------------------------------------------------------
    # Read colored boundary image
    # ---------------------------------------------------------

    boundary = cv2.imread(
        str(boundary_path),
        cv2.IMREAD_COLOR
    )

    if boundary is None:
        raise ValueError(
            f"Could not read boundary image:\n{boundary_path}"
        )

    # ---------------------------------------------------------
    # Make sure dimensions agree
    # ---------------------------------------------------------

    if filled.shape[:2] != boundary.shape[:2]:
        raise ValueError(
            f"Image dimensions do not match:\n"
            f"Filled:   {filled_path} -> {filled.shape[:2]}\n"
            f"Boundary: {boundary_path} -> {boundary.shape[:2]}"
        )

    # ---------------------------------------------------------
    # Find every non-black pixel in boundary image
    #
    # Because the boundaries may be red, yellow, blue, etc.,
    # check whether ANY of B, G, or R is > threshold.
    # ---------------------------------------------------------

    boundary_mask = np.any(
        boundary > boundary_threshold,
        axis=2
    )

    # ---------------------------------------------------------
    # Make filled mask strictly binary
    # ---------------------------------------------------------

    separated = (
        filled > 0
    ).astype(np.uint8) * 255

    # ---------------------------------------------------------
    # Remove boundary pixels
    # ---------------------------------------------------------

    separated[boundary_mask] = 0

    # ---------------------------------------------------------
    # Save result
    # ---------------------------------------------------------

    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    success = cv2.imwrite(
        str(output_path),
        separated
    )

    if not success:
        raise IOError(
            f"Could not save:\n{output_path}"
        )


def process_directory(
    filled_nuclei_dir,
    boundary_dir,
    output_dir,
    boundary_threshold=0,
):
    """
    Process all binary filled nuclei PNGs.

    Files are matched according to patch number.

    Example:

        epithelial_patch_7.png

    is matched with:

        overlay_patch_7.png
    """

    filled_nuclei_dir = Path(filled_nuclei_dir)
    boundary_dir = Path(boundary_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # ---------------------------------------------------------
    # Find all boundary images
    # ---------------------------------------------------------

    boundary_files = natsorted(
        boundary_dir.glob("*.png")
    )

    # Create dictionary:
    #
    # {
    #     1: overlay_patch_1.png,
    #     2: overlay_patch_2.png,
    #     ...
    # }
    boundary_lookup = {}

    for boundary_path in boundary_files:

        patch_number = extract_patch_number(
            boundary_path.name
        )

        if patch_number is not None:
            boundary_lookup[patch_number] = boundary_path

    print(
        f"Found {len(boundary_lookup)} boundary images."
    )

    # ---------------------------------------------------------
    # Find all filled nuclei images
    # ---------------------------------------------------------

    filled_files = natsorted(
        filled_nuclei_dir.glob("*.png")
    )

    print(
        f"Found {len(filled_files)} filled nuclei images."
    )

    processed = 0
    missing = 0

    # ---------------------------------------------------------
    # Process images
    # ---------------------------------------------------------

    for filled_path in filled_files:

        patch_number = extract_patch_number(
            filled_path.name
        )

        if patch_number is None:

            print(
                f"WARNING: Could not find patch number in "
                f"{filled_path.name}"
            )

            continue

        # Find corresponding boundary image
        boundary_path = boundary_lookup.get(
            patch_number
        )

        if boundary_path is None:

            print(
                f"WARNING: No boundary image found for "
                f"patch {patch_number}"
            )

            missing += 1
            continue

        # Keep original filled-mask filename
        output_path = (
            output_dir /
            filled_path.name
        )

        remove_boundaries(
            filled_path=filled_path,
            boundary_path=boundary_path,
            output_path=output_path,
            boundary_threshold=boundary_threshold,
        )

        processed += 1

        print(
            f"[{processed}/{len(filled_files)}] "
            f"{filled_path.name} "
            f"<- {boundary_path.name}"
        )

    print()
    print("Finished")
    print(f"Processed: {processed}")
    print(f"Missing boundary images: {missing}")


if __name__ == "__main__":

    # =========================================================
    # INPUT / OUTPUT PATHS
    # =========================================================

    # Example:
    # epithelial_patch_1.png
    # epithelial_patch_2.png
    # ...
    
    types = ["connective", "eosinophil", "epithelial", "lymphocyte", "neutrophil", "plasma"]
    
    for type in types:
        filled_nuclei_dir = (
            f"/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/with_only_binary_maps_from_saga/Func043_ST_HE_20x_BF_01_binary_only/overlay/binary_by_type/{type}/"
        )

        # Example:
        # overlay_patch_1.png
        # overlay_patch_2.png
        # ...
        boundary_dir = (
            "/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/Func043_ST_HE_20x_BF_01/overlay/"
        )

        # Processed masks will be saved here
        output_dir = (
            f"/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/Func043_ST_HE_20x_BF_01/binary_maps_separated_nuclei/{type}/"
        )

        # =========================================================
        # SETTINGS
        # =========================================================

        # Your uploaded boundary image has a perfectly black
        # background, so 0 is appropriate.
        #
        # Every pixel with ANY channel > 0 will be removed
        # from the filled mask.
        boundary_threshold = 0

        # =========================================================
        # RUN
        # =========================================================

        process_directory(
            filled_nuclei_dir=filled_nuclei_dir,
            boundary_dir=boundary_dir,
            output_dir=output_dir,
            boundary_threshold=boundary_threshold,
        )