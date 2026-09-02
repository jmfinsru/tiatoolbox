import pyvips
from PIL import Image
from pathlib import Path
import tifffile as tf

def convert_tiff(image_path : str | Path, save_path : str | Path):
    # Load the original TIFF file
    image = pyvips.Image.tiffload(image_path)
    # Save as a pyramidal tiled TIFF
    image.tiffsave(save_path,
                tile=True,
                pyramid=True,
                compression='deflate',  # or 'lzw', etc.
                tile_width=256,
                tile_height=256)  

def convert_and_resize_tiff(image_path : str | Path, save_path : str | Path):
    # Load the original TIFF file
    image = pyvips.Image.tiffload(image_path)
    # Resize image to 40x
    image = image.resize(2, kernel = "lanczos3") # Factor to scale image by and resampling kernel 
    # Save as a pyramidal tiled TIFF
    image.tiffsave(save_path,
                tile=True,
                pyramid=True,
                compression='deflate',  # or 'lzw', etc.
                tile_width=256,
                tile_height=256)  

from pathlib import Path
import pyvips

def convert_jpg_to_pyramid(image_path: str | Path, save_path: str | Path):
    image = pyvips.Image.new_from_file(str(image_path), access="sequential")

    image.tiffsave(
        str(save_path),
        tile=True,
        pyramid=True,
        subifd=True,        
        bigtiff=True,
        compression="deflate",
        level=6,
        tile_width=256,
        tile_height=256,
    )

    with tf.TiffFile(save_path) as tif:
        p0 = tif.pages[0]
        print("top-level pages:", len(tif.pages))
        print("page0 has SUBIFDs:", bool(p0.subifds))
        print("SUBIFD count (pyramid levels excluding base):", len(p0.subifds or ()))
    
    with tf.TiffFile(save_path) as tif:
        s0 = tif.series[0]
        print("levels:", len(s0.levels))
        for i, level in enumerate(s0.levels):
            page = level.pages[0]
            print(
                i,
                level.shape,
                "tiled:", page.is_tiled,
                "compression:", page.compression,
                "photometric:", page.photometric,
                "tile:", getattr(page, "tilewidth", None), getattr(page, "tilelength", None),
            )
if __name__ == "__main__":
    # Path to TIFF file 
    image_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF/HE_MIMB472_2B_120225.vsi.Collection/HE_MIMB472_2B_120225_20x_BF_01/HE_MIMB472_2B_120225_20x_BF_01.tif"
    save_path = "/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF_Pyramidal/Pyramidal_HE_MIMB472_2B_120225_20x_BF_01.tif"
    convert_tiff(image_path, save_path)
    # convert_jpg_to_pyramid(image_path, save_path)
    # convert_and_resize_tiff(image_path, save_path)
