from pathlib import Path
import pyvips


SUPPORTED_EXTENSIONS = {
    ".tif",
    ".tiff",
}


def make_output_path(
    image_path: Path,
    output_folder: Path,
) -> Path:
    output_name = image_path.name.replace(
        "_20x_",
        "_40x_upscale_",
    )

    return output_folder / output_name


def remove_metadata(image: pyvips.Image) -> pyvips.Image:
    """
    Remove metadata that can cause readers to interpret the image
    using stale OME-TIFF dimensions.
    """

    image = image.copy()

    metadata_fields_to_remove = [
        "image-description",   # often contains OME-XML
        "exif-data",
        "xmp-data",
        "iptc-data",
        "icc-profile-data",
        "photoshop-data",
    ]

    existing_fields = image.get_fields()

    for field in metadata_fields_to_remove:
        if field in existing_fields:
            image.remove(field)  # do not assign this back to image

    return image


def upscale_image(
    input_path: Path,
    output_path: Path,
    scale_factor: float = 2.0,
):
    image = pyvips.Image.new_from_file(
        str(input_path),
        access="sequential",
    )

    upscaled = image.resize(
        scale_factor,
        kernel="lanczos3",
    )

    upscaled = remove_metadata(upscaled)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # strip=True tells libvips not to copy metadata into the output TIFF
    upscaled.tiffsave(
        str(output_path),
        compression="lzw",
        tile=True,
        bigtiff=True,
        strip=True,
    )

    print(f"Saved: {output_path}")


def main(
    input_path: str,
    output_folder: str,
    scale_factor: float = 2.0,
):
    input_path = Path(input_path)
    output_folder = Path(output_folder)

    output_folder.mkdir(parents=True, exist_ok=True)

    if input_path.is_file():
        output_path = make_output_path(
            image_path=input_path,
            output_folder=output_folder,
        )

        upscale_image(
            input_path=input_path,
            output_path=output_path,
            scale_factor=scale_factor,
        )

    elif input_path.is_dir():
        image_paths = [
            path for path in input_path.iterdir()
            if path.suffix.lower() in SUPPORTED_EXTENSIONS
        ]

        if not image_paths:
            print(f"No supported TIFF images found in: {input_path}")
            return

        for image_path in sorted(image_paths):
            output_path = make_output_path(
                image_path=image_path,
                output_folder=output_folder,
            )

            upscale_image(
                input_path=image_path,
                output_path=output_path,
                scale_factor=scale_factor,
            )

    else:
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

if __name__ == "__main__":
    input_path = "/media/jenny/Expansion/HE/20x/images/Func116_ST_HE_20x_BF_01.tif"
    output_folder = "/media/jenny/Expansion/HE/40x_upscale/images/"
    scale_factor = 2.0

    main(
        input_path=input_path,
        output_folder=output_folder,
        scale_factor=scale_factor,
    )