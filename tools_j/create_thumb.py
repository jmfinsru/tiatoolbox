import pyvips
from pathlib import Path

def create_thumbnail(input_tiff_path : str | Path, output_thumbnail_path : str | Path, thumb_name: str):
    # Load the TIFF image
    image = pyvips.Image.new_from_file(input_tiff_path, access='sequential')

    # How much image is scaled
    scale = 0.2
    # Resize the image to the thumbnail size
    thumbnail = image.resize(scale)

    # Save thumbnail
    thumbnail.write_to_file(output_thumbnail_path / thumb_name)

if __name__ == "__main__":
    # Input and output path
    input_tiff_path = "/media/.../image.tif"
    output_thumbnail_path = Path("/media/.../")
    thumb_name = "x_scaled.png"
    if not output_thumbnail_path.exists():
        output_thumbnail_path.mkdir(parents=True) 
        print(f"Directory {output_thumbnail_path} was created")
    create_thumbnail(input_tiff_path, output_thumbnail_path, thumb_name)