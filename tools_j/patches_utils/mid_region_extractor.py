import glob
from natsort import natsorted
# from PIL import Image
# Image.MAX_IMAGE_PIXELS = None  # Set to None to disable the limit
from pathlib import Path
import pyvips
import matplotlib.pyplot as plt
from math import ceil
import numpy as np
from PIL import Image as pil_im

def create_wsi(tile_path: str | Path, patch_save_path: Path, wsi_path: str | Path):
    wsi_image = pyvips.Image.new_from_file(wsi_path)
    missing_pixels = 100
    # Regions in WSI that were cropped out when creating patches. Needs to be pasted into final image.
    region_wsi_lenght = wsi_image.crop(0, 0, wsi_image.width, missing_pixels)
    region_wsi_height = wsi_image.crop(0, 0, missing_pixels, wsi_image.height)

    # Find all images with .png
    image_list = glob.glob(tile_path + '*.tif')

    image_list = natsorted(image_list)
    # List with cropped regions of original image
    region_list = []

    # Loop through images
    for i in range(len(image_list)):
        image = pyvips.Image.new_from_file(image_list[i])
        # Define the region to extract (x, y, width, height)
        x = 100           # X coordinate of the left bound
        y = 100           # Y coordinate of the upper bound
        width = 1948     # Width of the region
        height = 1948    # Height of the region
        
        number_x = ceil(wsi_image.width/width)    # Number of patches along x axis
        number_y = ceil(wsi_image.height/height)   # Number of patches along y axis
        region = image.crop(x, y, width, height)
        region_list.append(region)
    
        ## Save cropped patches 
        # region.write_to_file("/media/jenny/Expansion/MM_HE_results/HE_HE_MM146_D_290125_20x_BF_01/wsi/tiles/" + f"extracted_region_{i}.png")

    # Create a new image with a size to fit all tiles
    new_image = pyvips.Image.black(region_list[0].width * number_x + missing_pixels, region_list[0].height * number_y + missing_pixels, bands = 3)
    # # Insert leftmost and uppermost cropped out pixels from original wsi file
    new_image = new_image.insert(region_wsi_lenght, 0, 0)
    new_image = new_image.insert(region_wsi_height, 0, 0)

    # Loop through the images and paste them into the new image
    for idx, img in enumerate(region_list):
        # Calculate the x and y position for the current image
        x_position = (idx % number_x) * region_list[0].width + missing_pixels
        y_position = (idx // number_x) * region_list[0].height + missing_pixels
        new_image=new_image.insert(img, x_position, y_position)
    
    remove_pixels_right = new_image.width - wsi_image.width # Pixels to remove from the right
    # remove_pixels_right = 0
    remove_pixels_bottom = new_image.height - wsi_image.height # Pixels to remove from the bottom
    # print(new_image.width)
    # print(new_image.height)
    # print(wsi_image.width)
    # print(wsi_image.height)
    # print(remove_pixels_bottom)
    # print(remove_pixels_right)

    # Calculate the new size after removing specified pixels from right and bottom
    new_image_width = new_image.width - remove_pixels_right
    new_image_height = new_image.height - remove_pixels_bottom
    print(f"Output image width: {new_image_width}")
    print(f"Output image height: {new_image_height}")
    # Crop the image
    new_image_cropped = new_image.crop(0, 0, new_image_width, new_image_height)

    # Save final result
    if not patch_save_path.exists():
        patch_save_path.mkdir(parents=True) 
        print(f"Directory {patch_save_path} was created")
    new_image_cropped.write_to_file( patch_save_path / "whole_image_complete.tif")



def create_binary_wsi(tile_path: str | Path, patch_save_path: Path, wsi_path: str | Path):
    tile_path = Path(tile_path)
    patch_save_path = Path(patch_save_path)
    patch_save_path.mkdir(parents=True, exist_ok=True)

    wsi_image = pyvips.Image.new_from_file(str(wsi_path), access="sequential")

    pixels = 3982
    nr_extra_pixels = 0
    patch_size = pixels + nr_extra_pixels  

    image_list = natsorted(
        list(tile_path.rglob("*Blood_vessel_binary_mask*.tif"))
    )

    print(f"Found {len(image_list)} mask patches")

    number_x = ceil(wsi_image.width / patch_size)
    number_y = ceil(wsi_image.height / patch_size)

    print("number_x:", number_x)
    print("number_y:", number_y)

    # Binary mask output: one band, black background
    new_image = pyvips.Image.black(
        wsi_image.width,
        wsi_image.height,
        bands=1
    ).cast("uchar")

    for idx, path in enumerate(image_list):
        tile = pyvips.Image.new_from_file(str(path), access="sequential")

        # Make one-band mask
        if tile.bands > 1:
            tile = tile.extract_band(0)

        # Convert to actual 0/255 values
        tile = (tile > 0).ifthenelse(255, 0).cast("uchar")

        col = idx % number_x
        row = idx // number_x

        x_position = col * patch_size
        y_position = row * patch_size

        # Crop tile at right/bottom WSI boundary
        insert_w = min(tile.width, wsi_image.width - x_position)
        insert_h = min(tile.height, wsi_image.height - y_position)

        if insert_w <= 0 or insert_h <= 0:
            continue

        tile = tile.crop(0, 0, insert_w, insert_h)

        new_image = new_image.insert(tile, x_position, y_position)

    out_path = patch_save_path / "Func040_VascMask.tif"

    new_image.write_to_file(
        str(out_path)
    )

    print("Saved:", out_path)

def find_png_files(collection_dir: str | Path) -> Path:
    collection_dir = Path(collection_dir)
    base_name = collection_dir.name
    collection_dir = collection_dir / "overlay/"
    print(f"collection_dir: {collection_dir}")
    

    # Find all.png files under this collection folder
    matches = natsorted(collection_dir.rglob(f"*.png"))

    if not matches:
        raise FileNotFoundError(f"No matching .png found in {collection_dir}")

    # # Prefer the standard brightfield image if it exists
    # preferred_name = f"{base_name}_20x_BF_0.tif"
    # preferred = [p for p in matches if p.name == preferred_name]
    # if preferred:
    #     return preferred[0]

    # Otherwise return the first matching tif
    return base_name, matches

def create_wsi_all_folders(root_dir: Path, root_wsi_dir: Path):

    # Find all images with .png
    for collection_dir in sorted(root_dir.glob("Func0*")):
        base_name , image_list = find_png_files(collection_dir)
    
        # Skip if this image has already been analyzed
        patch_save_path = Path(f"/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/{base_name}/wsi/")
        # # Skip if this image has already been analyzed
        # patch_save_path = Path(f"/media/jenny/Expansion/debugger/mid_region_extractor_debugger/{base_name}/wsi/")
    
        # if patch_save_path.exists():
        #     print(f"Skipping {base_name}, output already exists: {patch_save_path}")
        #     continue   # skip this image, go to next one
        if not patch_save_path.exists():
            patch_save_path.mkdir(parents=True)
            print(f"Directory {patch_save_path} was created")
        
        # Find WSI
        wsi_dir = root_wsi_dir
        # Find all tif files under this collection folder whose name starts with the base name
        wsi_path = next(wsi_dir.rglob(f"{base_name}_black_mask.tif"), None)
        
        if wsi_path is None:
            raise FileNotFoundError(f"No matching .tif found in {wsi_dir}")
        # Load WSI
        wsi_image = pyvips.Image.new_from_file(wsi_path)
        
        missing_pixels = 50
        # Regions in WSI that were cropped out when creating patches. Needs to be pasted into final image.
        region_wsi_lenght = wsi_image.crop(0, 0, wsi_image.width, missing_pixels)
        region_wsi_height = wsi_image.crop(0, 0, missing_pixels, wsi_image.height)
        # List with cropped regions of original image
        region_list = []
        
        # Patch index for pixel comparison between input patch and the same region in final WSI
        compare_idx = 2
        compare_region = None
        
        # Load and crop all input patches
        for i in range(len(image_list)):
            
            # Load image patch
            image = pyvips.Image.new_from_file(image_list[i])
        
            # Define the region to extract (x, y, width, height)
            x = 50           # X coordinate of the left bound
            y = 50           # Y coordinate of the upper bound
            width = 1948     # Width of the region
            height = 1948    # Height of the region
            
            number_x = ceil(wsi_image.width/width)    # Number of patches along x axis
            number_y = ceil(wsi_image.height/height)   # Number of patches along y axis
            region = image.crop(x, y, width, height)
            region_list.append(region)
            
            # Transform cropped region to a numpy array for comparison with final WSI
            if i == compare_idx:
                print(image_list[i])
                compare_region = region.numpy()
                if compare_region.ndim == 3 and compare_region.shape[2] == 4:
                    compare_region = compare_region[:, :, :3]

                
            ## Save cropped patches 
            # region.write_to_file("/media/jenny/Expansion/MM_HE_results/HE_HE_MM146_D_290125_20x_BF_01/wsi/tiles/" + f"extracted_region_{i}.png")
        # Create a new image with a size to fit all tiles
        new_image = pyvips.Image.black(region_list[0].width * number_x + missing_pixels, region_list[0].height * number_y + missing_pixels, bands = 3)
        
        # Insert leftmost and uppermost cropped out pixels from original wsi file
        new_image = new_image.insert(region_wsi_lenght, 0, 0)
        new_image = new_image.insert(region_wsi_height, 0, 0)

        # Loop through the image patches and paste them into the new image
        for idx, img in enumerate(region_list):
            if img.bands == 4:
                img = img.extract_band(0, n=3)
            # Calculate the x and y position for the current image
            x_position = (idx % number_x) * region_list[0].width + missing_pixels
            y_position = (idx // number_x) * region_list[0].height + missing_pixels
            new_image = new_image.insert(img, x_position, y_position)
        
        remove_pixels_right = new_image.width - wsi_image.width # Pixels to remove from the right
        # remove_pixels_right = 0
        remove_pixels_bottom = new_image.height - wsi_image.height # Pixels to remove from the bottom
        # print(new_image.width)
        # print(new_image.height)
        # print(wsi_image.width)
        # print(wsi_image.height)
        # print(remove_pixels_bottom)
        # print(remove_pixels_right)

        # Calculate the new size after removing specified pixels from right and bottom
        new_image_width = new_image.width - remove_pixels_right
        new_image_height = new_image.height - remove_pixels_bottom
        print(f"Output image width: {new_image_width}")
        print(f"Output image height: {new_image_height}")
        # Crop the image
        new_image_cropped = new_image.crop(0, 0, new_image_width, new_image_height)
        
        # Extract the exact location where tile idx=2 was inserted
        x0 = missing_pixels + (compare_idx % number_x) * width
        y0 = missing_pixels + (compare_idx // number_x) * height

        merged_region = new_image_cropped.crop(x0, y0, width, height).numpy()
        
        diff = merged_region.astype(np.int16) - compare_region.astype(np.int16)
        # print(compare_region[1500:1510, 1500:1510])
        # print("-----------------------------")
        # print(merged_region[1500:1510, 1500:1510])
        # print("-----------------------------")
        print("compare_region shape:", compare_region.shape)
        print("merged_region shape :", merged_region.shape)
        print("exact equal:", np.array_equal(merged_region, compare_region))
        print("diff min:", diff.min())
        print("diff max:", diff.max())
        print("max abs diff:", np.abs(diff).max())

        mask = np.any(diff != 0, axis=2)
        print("different pixels:", mask.sum())

        ys, xs = np.where(mask)
        for y, x in zip(ys[:10], xs[:10]):
            print(
                f"({y}, {x}) compare={compare_region[y, x]} "
                f"merged={merged_region[y, x]} diff={diff[y, x]}"
            )
        # Save final result
        if not patch_save_path.exists():
            patch_save_path.mkdir(parents=True) 
            print(f"Directory {patch_save_path} was created")
        new_image_cropped.write_to_file( patch_save_path / "whole_tp_image_complete.png")

        

if __name__ == "__main__":
    # # Path for output
    patch_save_path = Path("/media/jenny/Expansion/Prostata_QuPath/CD31_vessel_detection/Func040/binary_map/")
    # # Path to tiles
    tile_path = "/media/jenny/Expansion/Prostata_QuPath/CD31_vessel_detection/Func040/"
    # # Path to wsi
    wsi_path = "/media/jenny/Expansion/Prostata_Vilde/Analysis/Jenny/CD31_new_batch/10x/images/Func040_CD31.tif"
    # Path to root dir
    root_dir = Path("/media/jenny/Expansion/jenny_funcprost/conic/results/correct_nuclei/")
    # root_dir = Path("/media/jenny/Expansion/debugger/mid_region_extractor_debugger/")
    # Path to WSI root dir
    root_wsi_dir = Path("/media/jenny/Expansion/HE/20x/masks/")
    # Call function
    # create_wsi(tile_path, patch_save_path, wsi_path)
    # create_binary_wsi(tile_path, patch_save_path, wsi_path)
    create_wsi_all_folders(root_dir=root_dir, root_wsi_dir=root_wsi_dir)


