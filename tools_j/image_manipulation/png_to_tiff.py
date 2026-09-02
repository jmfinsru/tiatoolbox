import numpy as np
from PIL import Image
import tifffile as tiff
from pathlib import Path
import os


def png_to_pyramidal_ome_tiff(input_path, output_tiff_path, tile_size=512, compression="jpeg"):

            # Load PNG
            img = Image.open(input_path)
            img = img.convert("RGB")  
            base = np.array(img)
            
            # Create output path if it doesn't exist
            if not output_tiff_path.exists(): 
                output_tiff_path.mkdir(parents=True)
                print(f"Directory {output_tiff_path} was created")
            output_tiff = output_tiff_path
            print(output_tiff)

            # Build image pyramid
            pyramid = [base]
            while min(pyramid[-1].shape[0], pyramid[-1].shape[1]) > tile_size:
                down = Image.fromarray(pyramid[-1]).resize(
                    (pyramid[-1].shape[1] // 2, pyramid[-1].shape[0] // 2),
                    resample=Image.BILINEAR
                )
                pyramid.append(np.array(down))

            # Write pyramidal OME-TIFF
            with tiff.TiffWriter(output_tiff, bigtiff=True, ome=True) as tif:
                tif.write(
                    pyramid[0],
                    subifds=len(pyramid) - 1,
                    tile=(tile_size, tile_size),
                    compression=compression,
                    photometric="rgb",
                    metadata={
                        "axes": "YXS"
                    }
                )

                for level in pyramid[1:]:
                    tif.write(
                        level,
                        tile=(tile_size, tile_size),
                        compression=compression,
                        photometric="rgb"
                    )

            print(f"Saved pyramidal OME-TIFF: {output_tiff}")

    

if __name__ == "__main__":
    png_to_pyramidal_ome_tiff(
        input_path = "/media/jenny/Expansion/Arno/Image_M#29_CD4_Envision_tris-EDTA-rabbit.vsi.Collection/Image_M#29_CD4_Envision_tris-EDTA-rabbit_20x_BF_01/Image_M#29_CD4_Envision_tris-EDTA-rabbit_20x_BF_01.jpg",
        output_tiff_path = "/media/jenny/Expansion/Arno/Pyramidal/Pyramidal_Image_M#29_CD4_Envision_tris-EDTA-rabbit_20x_BF_01.ome.tif"
        )
