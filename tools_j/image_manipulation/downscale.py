import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

import tifffile
import pyvips


def get_pixels_element(root: ET.Element) -> ET.Element:
    """Find the first OME <Pixels> element."""
    for elem in root.iter():
        if elem.tag.endswith("Pixels"):
            return elem

    raise ValueError("No <Pixels> element found in OME-XML metadata")


def read_ome_metadata(tif_path: str | Path):
    """Read OME-XML and extract PhysicalSizeX/Y in µm/pixel."""
    tif_path = Path(tif_path)

    with tifffile.TiffFile(tif_path) as tif:
        omexml = tif.ome_metadata

    if omexml is None:
        raise ValueError(f"No OME-XML metadata found in: {tif_path}")

    root = ET.fromstring(omexml)

    # Preserve the default OME namespace when writing XML back
    if root.tag.startswith("{"):
        namespace = root.tag.split("}")[0].strip("{")
        ET.register_namespace("", namespace)

    pixels = get_pixels_element(root)

    physical_size_x = float(pixels.attrib["PhysicalSizeX"])
    physical_size_y = float(pixels.attrib["PhysicalSizeY"])

    return omexml, root, pixels, physical_size_x, physical_size_y


def update_ome_metadata_for_downscale(
    original_omexml: str,
    root: ET.Element,
    pixels: ET.Element,
    new_width: int,
    new_height: int,
    scale: float,
):
    """
    Update OME metadata after downscaling.

    If image is downscaled by 0.5:
        PhysicalSizeX: 0.25 -> 0.50 µm/pixel
        PhysicalSizeY: 0.25 -> 0.50 µm/pixel
    """

    old_physical_size_x = float(pixels.attrib["PhysicalSizeX"])
    old_physical_size_y = float(pixels.attrib["PhysicalSizeY"])

    new_physical_size_x = old_physical_size_x / scale
    new_physical_size_y = old_physical_size_y / scale

    pixels.attrib["SizeX"] = str(new_width)
    pixels.attrib["SizeY"] = str(new_height)

    pixels.attrib["PhysicalSizeX"] = f"{new_physical_size_x:.12g}"
    pixels.attrib["PhysicalSizeY"] = f"{new_physical_size_y:.12g}"

    # Keep units if already present; otherwise add micrometers
    pixels.attrib.setdefault("PhysicalSizeXUnit", "µm")
    pixels.attrib.setdefault("PhysicalSizeYUnit", "µm")

    updated_omexml = ET.tostring(root, encoding="unicode")

    # Add XML declaration back if the original had one
    if original_omexml.lstrip().startswith("<?xml"):
        updated_omexml = '<?xml version="1.0" encoding="UTF-8"?>\n' + updated_omexml

    return updated_omexml, new_physical_size_x, new_physical_size_y


def is_ome_tiff(path: str | Path) -> bool:
    name = Path(path).name.lower()
    return name.endswith(".ome.tif") or name.endswith(".ome.tiff")


def is_tiff(path: str | Path) -> bool:
    name = Path(path).name.lower()
    return (
        name.endswith(".tif")
        or name.endswith(".tiff")
        or name.endswith(".ome.tif")
        or name.endswith(".ome.tiff")
    )


def save_image(resized, output_path: Path):
    """
    Save image using the correct writer based on output extension.
    """

    suffix = output_path.suffix.lower()

    if is_tiff(output_path):
        resized.tiffsave(
            str(output_path),
            compression="lzw",
            tile=True,
            tile_width=512,
            tile_height=512,
            bigtiff=True
        )
        

    elif suffix in [".jpg", ".jpeg"]:
        resized.jpegsave(
            str(output_path),
            Q=95
        )
        print("Note: JPEG does not contain OME PhysicalSizeX/Y metadata.")
        print("Pixel-size metadata was not updated.")

    elif suffix == ".png":
        resized.pngsave(str(output_path))

    else:
        # Let pyvips infer from extension
        resized.write_to_file(str(output_path))


def downscale_image_half(
    input_path: str | Path,
    output_path: str | Path,
    is_mask: bool = False,
):
    input_path = Path(input_path)
    output_path = Path(output_path)

    scale = 0.5
    kernel = "nearest" if is_mask else "lanczos3"
    
    # Downscale image using pyvips
    img = pyvips.Image.new_from_file(str(input_path), access="sequential")
    resized = img.resize(scale, kernel=kernel)


    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    if input_path.suffix.lower() in [".tif", ".tiff", ".ome.tif", ".ome.tiff"]:

        # Read original OME metadata using tifffile
        original_omexml, root, pixels, old_psx, old_psy = read_ome_metadata(input_path)

        # Update OME metadata: dimensions and physical pixel size
        updated_omexml, new_psx, new_psy = update_ome_metadata_for_downscale(
            original_omexml=original_omexml,
            root=root,
            pixels=pixels,
            new_width=resized.width,
            new_height=resized.height,
            scale=scale,
        )

        # Attach updated OME-XML to TIFF ImageDescription
        resized = resized.copy()
        resized.set_type(
            pyvips.GValue.gstr_type,
            "image-description",
            updated_omexml
        )

        print(f"Input:  {input_path}")
        print(f"Output: {output_path}")
        print()
        print(f"Original image size: {img.width} x {img.height}")
        print(f"New image size:      {resized.width} x {resized.height}")
        print()
        print("Original pixel size:")
        print(f"  PhysicalSizeX = {old_psx} µm/pixel")
        print(f"  PhysicalSizeY = {old_psy} µm/pixel")
        print()
        print("New pixel size:")
        print(f"  PhysicalSizeX = {new_psx} µm/pixel")
        print(f"  PhysicalSizeY = {new_psy} µm/pixel")

    save_image(resized, output_path)


    
def verify_output_metadata(output_path: str | Path):
    """Optional check: read the output OME metadata and print PhysicalSizeX/Y."""
    _, _, pixels, psx, psy = read_ome_metadata(output_path)

    print()
    print("Verification from output OME metadata:")
    print(f"  SizeX = {pixels.attrib.get('SizeX')}")
    print(f"  SizeY = {pixels.attrib.get('SizeY')}")
    print(f"  PhysicalSizeX = {psx} µm/pixel")
    print(f"  PhysicalSizeY = {psy} µm/pixel")


def simple_downscale(input_path, output_path):

    scale = 0.75

    image = pyvips.Image.new_from_file(input_path, access="sequential")
    resized = image.resize(scale)

    resized.write_to_file(output_path)

    print(f"{image.width}x{image.height} -> {resized.width}x{resized.height}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Downsample an image by 0.5 while correcting pixel-size metadata."
    )
    parser.add_argument("--input", default="/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/Func050_ST_HE_20x_BF_01/wsi/whole_tp_image_complete.png", help="Input image path")
    parser.add_argument("--output", default="/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/Func050_ST_HE_20x_BF_01/wsi/whole_tp_image_complete_dowsampled_005.png", help="Output TIFF path")
    
    parser.add_argument(
        "--mask",
        action="store_true",
        help="Use nearest-neighbor interpolation for masks / label images"
    )

    args = parser.parse_args()

    # downscale_image_half(
    #     input_path=args.input,
    #     output_path=args.output,
    #     is_mask=args.mask,
    # )

    # verify_output_metadata(args.output)

    input_path = "/media/jenny/Expansion/HE/20x/images/Func116_ST_HE_20x_BF_01.tif"
    # output_path = "/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/Func050_ST_HE_20x_BF_01/wsi/connective/whole_connective_image_complete_downsampled_one_tenth.png"
    output_path = "/media/jenny/Expansion/HE/3_4_20x/images/Func116_ST_HE_20x_BF_01_dowsampled_three_fourth.png"
    simple_downscale(input_path, output_path)