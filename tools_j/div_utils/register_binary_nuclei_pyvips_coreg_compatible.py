from pathlib import Path
import csv
import gc

import numpy as np
import pyvips


# ============================================================
# MATRIX HELPERS
# ============================================================

def ensure_3x3_affine(matrix: np.ndarray) -> np.ndarray:
    """
    Convert 2x3 or 3x3 affine matrix to homogeneous 3x3 form.
    """
    matrix = np.asarray(matrix, dtype=np.float64)

    if matrix.shape == (2, 3):
        output = np.eye(3, dtype=np.float64)
        output[:2, :] = matrix
        matrix = output

    elif matrix.shape != (3, 3):
        raise ValueError(
            f"Expected a 2x3 or 3x3 affine matrix, got {matrix.shape}"
        )

    if not np.allclose(
        matrix[2, :],
        np.array([0.0, 0.0, 1.0]),
        atol=1e-10,
    ):
        raise ValueError(
            "Matrix is not an affine homogeneous transform.\n"
            f"Last row: {matrix[2, :]}"
        )

    return matrix


def load_final_matrix(
    final_matrix_path: str | Path,
) -> np.ndarray:
    """
    Load the final matrix produced by the co-registration script.

    IMPORTANT
    ---------
    The co-registration script defines:

        manual_matrix = rotation_translation @ scale_matrix

    and then:

        final_matrix = fine_residual_matrix @ manual_matrix

    Therefore final_matrix ALREADY contains:

        moving -> dimension scaling
               -> manual rotation/translation
               -> fine residual affine
               -> fixed H&E

    Do NOT multiply scale_matrix or manual_matrix into it again.
    """
    final_matrix = ensure_3x3_affine(
        np.load(final_matrix_path)
    )

    if not np.all(np.isfinite(final_matrix)):
        raise ValueError(
            "Final matrix contains NaN or infinite values."
        )

    determinant = np.linalg.det(
        final_matrix[:2, :2]
    )

    if abs(determinant) < 1e-12:
        raise ValueError(
            "Final registration matrix is singular."
        )

    print()
    print("=" * 70)
    print("FINAL CO-REGISTRATION MATRIX")
    print("=" * 70)
    print(final_matrix)
    print()
    print(f"Determinant: {determinant:.8f}")

    return final_matrix


def build_binary_to_coreg_moving_matrix(
    binary_width: int,
    binary_height: int,
    coreg_moving_width: int,
    coreg_moving_height: int,
) -> np.ndarray:
    """
    Convert coordinates from the full-resolution binary nuclei WSI
    into the coordinate system of the MOVING H&E image that was
    actually used by the co-registration script.

    Example:

        binary nuclei map:
            33654 x 32876

        co-registration moving H&E:
            approximately 3/4 of that size

    This conversion is separate from final_matrix.

    final_matrix expects coordinates belonging to the co-registration
    moving image, not the larger binary WSI.
    """

    scale_x = (
        coreg_moving_width
        /
        binary_width
    )

    scale_y = (
        coreg_moving_height
        /
        binary_height
    )

    matrix = np.array(
        [
            [scale_x, 0.0, 0.0],
            [0.0, scale_y, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )

    print()
    print("Binary -> co-registration moving-image matrix:")
    print(matrix)

    print(
        f"Binary -> moving scale X: {scale_x:.8f}"
    )
    print(
        f"Binary -> moving scale Y: {scale_y:.8f}"
    )

    return matrix


def build_binary_to_fixed_matrix(
    binary_width: int,
    binary_height: int,
    coreg_moving_width: int,
    coreg_moving_height: int,
    final_matrix: np.ndarray,
) -> np.ndarray:
    """
    Build the correct binary-WHI -> fixed-H&E transform.

    Coordinate chain:

        binary nuclei WSI
              |
              | binary_to_coreg_moving
              v
        exact moving-image coordinates used during co-registration
              |
              | final_matrix.npy
              v
        fixed H&E

    Therefore, using column-vector convention:

        binary_to_fixed =
            final_matrix @ binary_to_coreg_moving
    """

    binary_to_moving = (
        build_binary_to_coreg_moving_matrix(
            binary_width=binary_width,
            binary_height=binary_height,
            coreg_moving_width=coreg_moving_width,
            coreg_moving_height=coreg_moving_height,
        )
    )

    binary_to_fixed = (
        ensure_3x3_affine(final_matrix)
        @
        binary_to_moving
    )

    print()
    print("Complete BINARY -> FIXED matrix:")
    print(binary_to_fixed)

    a = binary_to_fixed[0, 0]
    b = binary_to_fixed[0, 1]
    c = binary_to_fixed[1, 0]
    d = binary_to_fixed[1, 1]

    effective_scale_x = np.sqrt(
        a ** 2 + c ** 2
    )

    effective_scale_y = np.sqrt(
        b ** 2 + d ** 2
    )

    print()
    print(
        f"Effective binary -> fixed scale X: "
        f"{effective_scale_x:.8f}"
    )
    print(
        f"Effective binary -> fixed scale Y: "
        f"{effective_scale_y:.8f}"
    )
    print(
        f"Effective determinant: "
        f"{np.linalg.det(binary_to_fixed[:2, :2]):.8f}"
    )

    return binary_to_fixed


# ============================================================
# IMAGE HELPERS
# ============================================================

def get_image_dimensions(
    image_path: str | Path,
) -> tuple[int, int]:
    """
    Return width, height without converting a WSI to NumPy.
    """
    image = pyvips.Image.new_from_file(
        str(image_path),
        access="sequential",
    )

    width = image.width
    height = image.height

    del image

    return width, height


def load_binary_vips(
    image_path: str | Path,
    access: str = "random",
) -> pyvips.Image:
    """
    Load binary nuclei map and force strict one-band uint8 0/255.
    """
    image = pyvips.Image.new_from_file(
        str(image_path),
        access=access,
    )

    if image.bands > 1:
        image = image.extract_band(0)

    image = (
        image > 0
    ).ifthenelse(
        255,
        0,
    ).cast(
        "uchar"
    )

    return image


def load_fixed_rgb_vips(
    fixed_image_path: str | Path,
) -> pyvips.Image:
    """
    Load fixed H&E as 3-band uchar.
    """
    image = pyvips.Image.new_from_file(
        str(fixed_image_path),
        access="random",
    )

    if image.bands >= 3:
        image = image.extract_band(
            0,
            n=3,
        )

    elif image.bands == 1:
        image = image.bandjoin(
            [image, image]
        )

    elif image.bands == 2:
        first = image.extract_band(0)
        image = first.bandjoin(
            [first, first]
        )

    if image.format != "uchar":
        image = image.cast("uchar")

    return image


def save_vips_image(
    image: pyvips.Image,
    output_path: str | Path,
) -> None:
    """
    Save WSI-sized image without conversion to NumPy.
    """
    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    suffix = output_path.suffix.lower()

    if suffix == ".png":

        image.pngsave(
            str(output_path),
            compression=6,
        )

    elif suffix in {
        ".tif",
        ".tiff",
    }:

        image.tiffsave(
            str(output_path),
            compression="deflate",
            tile=True,
            tile_width=512,
            tile_height=512,
            bigtiff=True,
            pyramid=False,
        )

    else:

        image.write_to_file(
            str(output_path)
        )

    print(
        f"Saved:\n{output_path}"
    )


# ============================================================
# DIAGNOSTICS
# ============================================================

def inspect_transformed_bounds(
    matrix: np.ndarray,
    source_width: int,
    source_height: int,
    fixed_width: int,
    fixed_height: int,
) -> None:
    """
    Print where the four binary-map corners land in fixed coordinates.
    """

    matrix = ensure_3x3_affine(
        matrix
    )

    corners = np.array(
        [
            [0, 0, 1],
            [source_width - 1, 0, 1],
            [0, source_height - 1, 1],
            [
                source_width - 1,
                source_height - 1,
                1,
            ],
        ],
        dtype=np.float64,
    ).T

    transformed = (
        matrix
        @
        corners
    )

    transformed = transformed[:2, :].T

    names = [
        "top-left",
        "top-right",
        "bottom-left",
        "bottom-right",
    ]

    print()
    print("=" * 70)
    print("TRANSFORMED BINARY-MAP CORNERS")
    print("=" * 70)

    for name, point in zip(
        names,
        transformed,
    ):
        print(
            f"{name:<12} "
            f"x={point[0]:10.2f}, "
            f"y={point[1]:10.2f}"
        )

    min_x = transformed[:, 0].min()
    max_x = transformed[:, 0].max()

    min_y = transformed[:, 1].min()
    max_y = transformed[:, 1].max()

    print()
    print(
        f"Transformed bounds:"
    )
    print(
        f"X: {min_x:.2f} -> {max_x:.2f}"
    )
    print(
        f"Y: {min_y:.2f} -> {max_y:.2f}"
    )

    print()
    print(
        f"Fixed canvas:"
    )
    print(
        f"X: 0 -> {fixed_width - 1}"
    )
    print(
        f"Y: 0 -> {fixed_height - 1}"
    )


# ============================================================
# FIND BINARY MAP
# ============================================================

def find_binary_type_map(
    type_dir: str | Path,
) -> Path:
    """
    Find one reconstructed whole-slide binary map for a nucleus type.

    Preferred:
        whole_<type>_image_complete.png
        whole_<type>_image_complete.tif
    """
    type_dir = Path(type_dir)
    type_name = type_dir.name

    preferred = [
        type_dir
        /
        f"whole_{type_name}_image_complete.png",

        type_dir
        /
        f"whole_{type_name}_image_complete.tif",

        type_dir
        /
        f"whole_{type_name}_image_complete.tiff",
    ]

    for path in preferred:

        if path.exists():
            return path

    matches = []

    for pattern in (
        "*.png",
        "*.tif",
        "*.tiff",
    ):
        matches.extend(
            type_dir.rglob(pattern)
        )

    matches = sorted(
        set(matches),
        key=lambda p: str(p),
    )

    if not matches:
        raise FileNotFoundError(
            f"No binary WSI found inside:\n"
            f"{type_dir}"
        )

    if len(matches) == 1:
        return matches[0]

    raise RuntimeError(
        f"Multiple images found in {type_dir}, "
        "but the expected whole-slide filename is absent."
    )


# ============================================================
# LIBVIPS AFFINE REGISTRATION
# ============================================================

def register_binary_map(
    binary_path: str | Path,
    registration_matrix: np.ndarray,
    fixed_width: int,
    fixed_height: int,
) -> pyvips.Image:
    """
    Transform a full-resolution binary nuclei WSI to fixed H&E
    coordinates with nearest-neighbour interpolation.
    """

    binary = load_binary_vips(
        binary_path,
        access="random",
    )

    matrix = ensure_3x3_affine(
        registration_matrix
    )

    a = float(matrix[0, 0])
    b = float(matrix[0, 1])
    tx = float(matrix[0, 2])

    c = float(matrix[1, 0])
    d = float(matrix[1, 1])
    ty = float(matrix[1, 2])

    print()
    print("Applying binary -> fixed affine:")
    print(
        f"[{a:.10f}  {b:.10f}  {tx:.10f}]"
    )
    print(
        f"[{c:.10f}  {d:.10f}  {ty:.10f}]"
    )

    nearest = pyvips.Interpolate.new(
        "nearest"
    )

    registered = binary.affine(
        [
            a,
            b,
            c,
            d,
        ],

        interpolate=nearest,

        odx=tx,
        ody=ty,

        # Exactly the fixed-H&E coordinate canvas.
        oarea=[
            0,
            0,
            int(fixed_width),
            int(fixed_height),
        ],

        background=[0],
        extend="background",
    )

    registered = (
        registered > 0
    ).ifthenelse(
        255,
        0,
    ).cast(
        "uchar"
    )

    return registered


# ============================================================
# OPTIONAL BINARY OVERLAY
# ============================================================

def create_binary_overlay(
    fixed_rgb: pyvips.Image,
    registered_binary: pyvips.Image,
    alpha: float = 0.50,
) -> pyvips.Image:
    """
    Overlay white nuclei on fixed H&E.
    """

    if not 0.0 <= alpha <= 1.0:
        raise ValueError(
            "alpha must be between 0 and 1"
        )

    foreground = (
        registered_binary > 0
    )

    whitened = fixed_rgb.linear(
        1.0 - alpha,
        255.0 * alpha,
    ).cast(
        "uchar"
    )

    overlay = foreground.ifthenelse(
        whitened,
        fixed_rgb,
    )

    return overlay


# ============================================================
# OPTIONAL PURE-LIBVIPS REGION COUNT
# ============================================================

def count_instances_vips(
    binary_map: pyvips.Image,
) -> dict:
    """
    Count 4-connected regions without converting the full WSI to NumPy.

    Assumes one connected black background region.

    This does not perform a minimum-area filter.
    """

    binary_map = (
        binary_map > 0
    ).ifthenelse(
        255,
        0,
    ).cast(
        "uchar"
    )

    label_result = binary_map.labelregions(
        segments=True
    )

    if not isinstance(
        label_result,
        list,
    ):
        raise RuntimeError(
            "Unexpected pyvips labelregions return value."
        )

    labelled = label_result[0]
    extra = label_result[1]

    if isinstance(extra, dict):
        total_regions = int(
            extra["segments"]
        )
    else:
        total_regions = int(
            extra
        )

    instance_count = max(
        total_regions - 1,
        0,
    )

    foreground_pixels = int(
        round(
            (
                float(binary_map.avg())
                /
                255.0
            )
            *
            binary_map.width
            *
            binary_map.height
        )
    )

    mean_area = (
        foreground_pixels
        /
        instance_count
        if instance_count > 0
        else 0.0
    )

    del labelled

    return {
        "instance_count":
            instance_count,

        "foreground_pixels":
            foreground_pixels,

        "mean_instance_area_px":
            float(mean_area),
    }


def save_instance_csv(
    results: list[dict],
    csv_path: str | Path,
) -> None:
    """
    Save optional instance statistics.
    """

    csv_path = Path(csv_path)

    csv_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = [
        "nuclei_type",
        "instance_count",
        "foreground_pixels",
        "mean_instance_area_px",
    ]

    with open(
        csv_path,
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for row in results:
            writer.writerow(row)


# ============================================================
# MAIN REGISTRATION FUNCTION
# ============================================================

def register_all_nuclei_maps(
    type_root: str | Path,
    fixed_image_path: str | Path,
    coreg_moving_image_path: str | Path,
    final_matrix_path: str | Path,
    save_dir: str | Path,
    output_extension: str = ".png",
    save_binary_overlay: bool = False,
    overlay_alpha: float = 0.50,
    count_instances: bool = False,
) -> None:
    """
    Register all binary nucleus-type maps in a structure such as:

        type_root/
            connective/
                whole_connective_image_complete.png
            epithelial/
                whole_epithelial_image_complete.png
            eosinophil/
                whole_eosinophil_image_complete.png
            lymphocyte/
                whole_lymphocyte_image_complete.png
            neutrophil/
                whole_neutrophil_image_complete.png
            plasma/
                whole_plasma_image_complete.png

    CRITICAL:
        coreg_moving_image_path must point to the EXACT moving H&E
        image used when final_matrix.npy was generated.

    The binary maps may be larger than that moving H&E. Their coordinates
    are first converted to the co-registration moving-image coordinate
    system, then final_matrix is applied.
    """

    type_root = Path(
        type_root
    )

    fixed_image_path = Path(
        fixed_image_path
    )

    coreg_moving_image_path = Path(
        coreg_moving_image_path
    )

    final_matrix_path = Path(
        final_matrix_path
    )

    save_dir = Path(
        save_dir
    )

    save_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not output_extension.startswith("."):
        output_extension = (
            "."
            +
            output_extension
        )

    # ========================================================
    # FIXED IMAGE
    # ========================================================

    fixed_width, fixed_height = (
        get_image_dimensions(
            fixed_image_path
        )
    )

    print()
    print("=" * 70)
    print("FIXED H&E")
    print("=" * 70)
    print(
        f"{fixed_width:,} x "
        f"{fixed_height:,}"
    )

    # ========================================================
    # EXACT MOVING IMAGE USED BY CO-REGISTRATION
    # ========================================================

    moving_width, moving_height = (
        get_image_dimensions(
            coreg_moving_image_path
        )
    )

    print()
    print("=" * 70)
    print("CO-REGISTRATION MOVING H&E")
    print("=" * 70)
    print(
        f"{moving_width:,} x "
        f"{moving_height:,}"
    )
    print(
        f"Path:\n{coreg_moving_image_path}"
    )

    # ========================================================
    # FINAL MATRIX
    # ========================================================

    final_matrix = (
        load_final_matrix(
            final_matrix_path
        )
    )

    # ========================================================
    # FIND TYPE MAPS
    # ========================================================

    type_dirs = sorted(
        [
            path
            for path in type_root.iterdir()
            if path.is_dir()
        ],
        key=lambda p: p.name,
    )

    if not type_dirs:
        raise FileNotFoundError(
            f"No nucleus-type directories inside:\n"
            f"{type_root}"
        )

    valid_maps = []

    for type_dir in type_dirs:

        try:
            binary_path = (
                find_binary_type_map(
                    type_dir
                )
            )

        except (
            FileNotFoundError,
            RuntimeError,
        ) as error:

            print(
                f"Skipping {type_dir.name}: "
                f"{error}"
            )

            continue

        valid_maps.append(
            (
                type_dir.name,
                binary_path,
            )
        )

    if not valid_maps:
        raise FileNotFoundError(
            "No usable binary nucleus maps found."
        )

    # ========================================================
    # OPTIONAL FIXED IMAGE
    # ========================================================

    fixed_rgb = None

    if save_binary_overlay:

        fixed_rgb = (
            load_fixed_rgb_vips(
                fixed_image_path
            )
        )

    # ========================================================
    # PROCESS TYPES
    # ========================================================

    count_results = []

    reference_binary_size = None

    for index, (
        type_name,
        binary_path,
    ) in enumerate(
        valid_maps,
        start=1,
    ):

        print()
        print("=" * 70)
        print(
            f"[{index}/{len(valid_maps)}] "
            f"PROCESSING: {type_name}"
        )
        print("=" * 70)

        binary_width, binary_height = (
            get_image_dimensions(
                binary_path
            )
        )

        print()
        print(
            f"Binary dimensions: "
            f"{binary_width:,} x "
            f"{binary_height:,}"
        )

        # All type maps should represent the same original WSI.
        if reference_binary_size is None:

            reference_binary_size = (
                binary_width,
                binary_height,
            )

        elif (
            binary_width,
            binary_height,
        ) != reference_binary_size:

            raise ValueError(
                "Nucleus-type binary maps do not all have the "
                "same dimensions.\n"
                f"Expected: {reference_binary_size}\n"
                f"{type_name}: "
                f"{binary_width} x {binary_height}"
            )

        # ----------------------------------------------------
        # THIS IS THE CRITICAL CORRECT COMPOSITION
        #
        # binary coordinates
        #   -> exact co-reg moving-image coordinates
        #   -> final_matrix
        #   -> fixed
        #
        # final_matrix itself ALREADY contains the co-reg
        # scale + manual transform + fine residual.
        # ----------------------------------------------------

        registration_matrix = (
            build_binary_to_fixed_matrix(
                binary_width=binary_width,
                binary_height=binary_height,
                coreg_moving_width=moving_width,
                coreg_moving_height=moving_height,
                final_matrix=final_matrix,
            )
        )

        inspect_transformed_bounds(
            matrix=registration_matrix,
            source_width=binary_width,
            source_height=binary_height,
            fixed_width=fixed_width,
            fixed_height=fixed_height,
        )

        registered_binary = (
            register_binary_map(
                binary_path=binary_path,
                registration_matrix=registration_matrix,
                fixed_width=fixed_width,
                fixed_height=fixed_height,
            )
        )

        output_path = (
            save_dir
            /
            f"{type_name}_registered_binary"
            f"{output_extension}"
        )

        save_vips_image(
            registered_binary,
            output_path,
        )

        del registered_binary
        gc.collect()

        # ----------------------------------------------------
        # Reopen written result only if another operation needs it
        # ----------------------------------------------------

        registered_saved = None

        if (
            save_binary_overlay
            or
            count_instances
        ):

            registered_saved = (
                load_binary_vips(
                    output_path,
                    access="random",
                )
            )

        if save_binary_overlay:

            overlay = (
                create_binary_overlay(
                    fixed_rgb=fixed_rgb,
                    registered_binary=registered_saved,
                    alpha=overlay_alpha,
                )
            )

            overlay_path = (
                save_dir
                /
                f"{type_name}_registered_binary_overlay"
                f"{output_extension}"
            )

            save_vips_image(
                overlay,
                overlay_path,
            )

            del overlay

        if count_instances:

            stats = (
                count_instances_vips(
                    registered_saved
                )
            )

            result = {
                "nuclei_type":
                    type_name,

                **stats,
            }

            count_results.append(
                result
            )

        if registered_saved is not None:
            del registered_saved

        gc.collect()

    # ========================================================
    # OPTIONAL CSV
    # ========================================================

    if count_instances:

        csv_path = (
            save_dir
            /
            "nuclei_instance_counts.csv"
        )

        save_instance_csv(
            count_results,
            csv_path,
        )

    if fixed_rgb is not None:
        del fixed_rgb

    gc.collect()

    print()
    print("=" * 70)
    print("FINISHED")
    print("=" * 70)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    biopsy_number = "116"

    # ========================================================
    # FIXED H&E
    # ========================================================

    fixed_image_path = Path(
        "/media/jenny/Expansion1/Prostata_Vilde/"
        "Co-reg_images_10x/"
        f"Func{biopsy_number}_HE.tif"
    )

    # ========================================================
    # EXACT MOVING H&E USED BY THE CO-REGISTRATION SCRIPT
    #
    # This must be the SAME file as moving_path in the
    # co-registration script that generated final_matrix.npy.
    #
    # For Func050 that script used:
    #
    # Func050_ST_HE_20x_BF_01_dowsampled_three_fourth.png
    # ========================================================

    coreg_moving_image_path = Path(
        "/media/jenny/Expansion/HE/"
        "3_4_20x/images/"
        f"Func{biopsy_number}_ST_HE_20x_BF_01_"
        "dowsampled_three_fourth.png"
    )

    # ========================================================
    # BINARY NUCLEUS-TYPE WHOLE-SLIDE MAPS
    # ========================================================

    type_root = Path(
        "/media/jenny/Expansion/jenny_funcprost/"
        "conic/results/correct_nuclei/"
        f"Func{biopsy_number}_ST_HE_20x_BF_01/"
        "wsi/"
    )

    # ========================================================
    # FINAL MATRIX PRODUCED BY CO-REGISTRATION
    #
    # Do NOT additionally multiply scale_matrix.npy or
    # manual_matrix.npy into this.
    #
    # The co-registration code has already composed them into
    # final_matrix.npy.
    # ========================================================

    final_matrix_path = Path(
        "/media/jenny/Expansion/jenny_funcprost/"
        "conic/results/correct_nuclei/"
        f"Func{biopsy_number}_ST_HE_20x_BF_01/"
        "co_registration/"
        "final_matrix.npy"
    )

    # ========================================================
    # OUTPUT
    # ========================================================

    save_dir = Path(
        "/media/jenny/Expansion/jenny_funcprost/"
        "conic/results/correct_nuclei/"
        f"Func{biopsy_number}_ST_HE_20x_BF_01/"
        "co_registration/"
        "registered_nuclei_maps"
    )

    # ========================================================
    # OPTIONS
    # ========================================================

    # For especially large output maps, ".tif" will write tiled BigTIFF.
    output_extension = ".png"

    save_binary_overlay = False

    overlay_alpha = 0.50

    # Region counting can be expensive for whole-slide maps.
    count_instances = False

    # ========================================================
    # RUN
    # ========================================================

    register_all_nuclei_maps(
        type_root=type_root,
        fixed_image_path=fixed_image_path,
        coreg_moving_image_path=coreg_moving_image_path,
        final_matrix_path=final_matrix_path,
        save_dir=save_dir,
        output_extension=output_extension,
        save_binary_overlay=save_binary_overlay,
        overlay_alpha=overlay_alpha,
        count_instances=count_instances,
    )
