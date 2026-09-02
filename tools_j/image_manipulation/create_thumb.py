import pyvips
from pathlib import Path

def create_thumbnail(input_tiff_path : str | Path, output_thumbnail_path : str | Path, thumb_name: str):
    # Load the TIFF image
    image = pyvips.Image.new_from_file(input_tiff_path, access='sequential')

    # How much image is scaled
    scale = 0.08
    # Resize the image to the thumbnail size
    thumbnail = image.resize(scale)

    # Save thumbnail
    thumbnail.write_to_file(output_thumbnail_path / thumb_name)

from pathlib import Path
import pyvips


def create_thumbnail_new(input_image_path: str | Path, output_thumbnail_path: str | Path, thumb_name: str, scale: float = 0.08):
    input_image_path = Path(input_image_path)
    output_thumbnail_path = Path(output_thumbnail_path)

    if not output_thumbnail_path.exists():
        output_thumbnail_path.mkdir(parents=True, exist_ok=True)

    # Load image
    image = pyvips.Image.new_from_file(str(input_image_path), access="sequential")

    # Resize
    thumbnail = image.resize(scale)

    # Save thumbnail
    thumbnail.write_to_file(str(output_thumbnail_path / thumb_name))


def create_thumbnails_all(root_dir: str | Path, output_thumbnail_path: str | Path, scale: float = 0.03):
    root_dir = Path(root_dir)
    output_thumbnail_path = Path(output_thumbnail_path)

    # Find all files named whole_image_complete.png inside wsi folders
    image_paths = sorted(root_dir.rglob("wsi/whole_image_complete.png"))

    if not image_paths:
        raise FileNotFoundError(f"No whole_image_complete.png files found under {root_dir}")

    print(f"Found {len(image_paths)} images")

    for image_path in image_paths:
        # Example:
        # /media/jenny/Expansion/HE_20x_results/HE_SIDB472_2B_120225/wsi/whole_image_complete.png
        # sample folder = HE_SIDB472_2B_120225
        sample_name = image_path.parent.parent.name

        thumb_name = f"{sample_name}_thumbnail.png"

        print(f"Creating thumbnail for: {image_path}")
        create_thumbnail_new(
            input_image_path=image_path,
            output_thumbnail_path=output_thumbnail_path,
            thumb_name=thumb_name,
            scale=scale,
        )

    print(f"Saved thumbnails in: {output_thumbnail_path}")


if __name__ == "__main__":
    #Input and output path
    input_tiff_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF/HE_MIMB711_2A_120225.vsi.Collection/HE_MIMB711_2A_120225_20x_BF_01/HE_MIMB711_2A_120225_20x_BF_01.tif"
    output_thumbnail_path = Path("/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF_thumb/")
    thumb_name = "HE_MIMB711_2A_120225_20x_BF_01_thumb.png"
    if not output_thumbnail_path.exists():
        output_thumbnail_path.mkdir(parents=True) 
        print(f"Directory {output_thumbnail_path} was created")
    create_thumbnail(input_tiff_path, output_thumbnail_path, thumb_name)

   
    # root_dir = Path("/media/jenny/Expansion/MetoxyLacc_HE_20x_results")
    # output_thumbnail_path = Path("/media/jenny/Expansion/MetoxyLacc_HE_20x_results/thumbnails")

    # create_thumbnails_all(root_dir, output_thumbnail_path, scale=0.08)
