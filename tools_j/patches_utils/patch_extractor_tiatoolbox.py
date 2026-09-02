import logging

if logging.getLogger().hasHandlers():
    logging.getLogger().handlers.clear()

import sys
import tifffile
import os
import numpy as np
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

# Get the path to the project root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Append the project root to sys.path
if project_root not in sys.path:
    sys.path.append(project_root)

from tiatoolbox.wsicore import WSIReader
from tiatoolbox.tools.patchextraction import SlidingWindowPatchExtractor


def find_main_tif(collection_dir: str | Path) -> Path:
    collection_dir = Path(collection_dir)
    print(f"collection_dir: {collection_dir}")
    base_name = collection_dir.name.removesuffix(".vsi.Collection")

    # Find all tif files under this collection folder whose name starts with the base name
    matches = sorted(collection_dir.rglob(f"{base_name}_*.tif"))

    if not matches:
        raise FileNotFoundError(f"No matching .tif found in {collection_dir}")

    # Prefer the standard brightfield image if it exists
    preferred_name = f"{base_name}_20x_BF_0.tif"
    preferred = [p for p in matches if p.name == preferred_name]
    if preferred:
        return preferred[0]

    # Otherwise return the first matching tif
    return base_name, matches[0]


def extract_patches_all_in_folder(root_dir : str | Path):

    for collection_dir in sorted(root_dir.glob("*.vsi.Collection")):
        base_name , image_path = find_main_tif(collection_dir)
        # Skip if this image has already been analyzed

        patch_save_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_patches/{base_name}/2048x2048/")
    
        if patch_save_path.exists():
            print(f"Skipping {base_name}, output already exists: {patch_save_path}")
            continue   # skip this image, go to next one
        
        
        patch_save_path.mkdir(parents=True)
        print(f"Directory {patch_save_path} was created")

        print(f"Processing: {image_path}")

        try:
            img = tifffile.imread(image_path)
        except:
            img = image_path
        
        nr_extra_pixels = 0        # Number of extra pixels per patch
        pixels = 2048               # Pixels along x- and y- axis for each patch before extra pixels are added. x and y must have same number of pixels if processed with aug_hovernet.
        nr_pixels_overlap = 100     # Number of pixels overlap between patches
        # number of pixels along x- and y axis for each patch
        x = pixels + nr_extra_pixels
        y = pixels + nr_extra_pixels

        
        extractor = SlidingWindowPatchExtractor(
            img,
            patch_size=(x, y),
            stride=(x - nr_pixels_overlap, y - nr_pixels_overlap),
            pad_mode="constant",
            pad_constant_values=255,   # white padding
            within_bound=False
        )
        
        print(f"Number of images: {len(extractor)}")

        # Iterate over the patches and save each patch
        for patch_idx in range(len(extractor)): 
            patch = extractor[patch_idx]
            patch_image = Image.fromarray(patch)
            patch_filename = os.path.join(patch_save_path, f"patch_{patch_idx}.png")
            patch_image.save(patch_filename)

def find_main_tif_mask(collection_dir: str | Path) -> Path:
    collection_dir = Path(collection_dir)
    print(f"collection_dir: {collection_dir}")
    base_name = collection_dir.name

    # Find specific tif files under this collection folder
    matches = sorted(collection_dir.rglob(f"{base_name}_no_folds_mask.tif"))

    if not matches:
        raise FileNotFoundError(f"No matching .tif found in {collection_dir}")

    # # Prefer the standard brightfield image if it exists
    # preferred_name = f"{base_name}_20x_BF_0.tif"
    # preferred = [p for p in matches if p.name == preferred_name]
    # if preferred:
    #     return preferred[0]

    # Otherwise return the first matching tif
    return base_name, matches[0]


def extract_patches_all_in_folder_mask(root_dir : str | Path):
    
    for collection_dir in sorted(root_dir.glob("HE_MIM*")):
        base_name , image_path = find_main_tif_mask(collection_dir)
        # Skip if this image has already been analyzed
        base_name = base_name.removesuffix("_20x_BF_01")
        patch_save_path = Path(f"/media/jenny/Expansion/MetoxyLacc_HE_20x_patches/{base_name}/2048x2048_mask/")
    
        if patch_save_path.exists():
            print(f"Skipping {base_name}, output already exists: {patch_save_path}")
            continue   # skip this image, go to next one
        
        
        patch_save_path.mkdir(parents=True)
        print(f"Directory {patch_save_path} was created")

        print(f"Processing: {image_path}")

        try:
            img = tifffile.imread(image_path)
        except:
            img = image_path
        
        nr_extra_pixels = 0        # Number of extra pixels per patch
        pixels = 2048               # Pixels along x- and y- axis for each patch before extra pixels are added. x and y must have same number of pixels if processed with aug_hovernet.
        nr_pixels_overlap = 100     # Number of pixels overlap between patches
        # number of pixels along x- and y axis for each patch
        x = pixels + nr_extra_pixels
        y = pixels + nr_extra_pixels

        
        extractor = SlidingWindowPatchExtractor(
            img,
            patch_size=(x, y),
            stride=(x - nr_pixels_overlap, y - nr_pixels_overlap),
            pad_mode="constant",
            pad_constant_values=255,   # white padding
            within_bound=False
        )
        
        print(f"Number of images: {len(extractor)}")

        # Iterate over the patches and save each patch
        for patch_idx in range(len(extractor)): 
            patch = extractor[patch_idx]
            patch_image = Image.fromarray(patch)
            patch_filename = os.path.join(patch_save_path, f"patch_{patch_idx}.png")
            patch_image.save(patch_filename)


def extract_patches_wsi(wsi_path : str | Path, patch_save_path : Path, mask_path: str | Path = None):


    wsi = WSIReader.open(wsi_path)
    dim = wsi.slide_dimensions(resolution=1, units="power")
    # Print pixels along x- and y axis for wsi
    print(f"wsi dim: {dim}")

    nr_extra_pixels = 0          # Number of extra pixels per patch
    pixels = 2048                # Pixels along x- and y- axis for each patch before extra pixels are added. x and y must have same number of pixels in numbers of 2^x if processed with aug_hovernet.
    nr_pixels_overlap = 100     # Number of pixels overlap between patches
    # number of pixels along x- and y axis for each patch
    x = pixels + nr_extra_pixels
    y = pixels + nr_extra_pixels

    extractor = SlidingWindowPatchExtractor(
        wsi_path,
        patch_size=(x, y),
        stride=(x-nr_pixels_overlap, y-nr_pixels_overlap),
        # mask_path

    )
    print(f"Number of images: {len(extractor)}")

    # Iterate over the patches and save each patch
    for patch_idx in range(len(extractor)): 
        patch = extractor[patch_idx]

        patch_image = Image.fromarray(patch)
        patch_filename = os.path.join(patch_save_path, f"patch_{patch_idx}.png")
        patch_image.save(patch_filename)

def extract_patches(image_path : str | Path, patch_save_path : Path):
    
    try:
        img = tifffile.imread(image_path)
    except:
        img = image_path
    
    nr_extra_pixels = 0        # Number of extra pixels per patch
    pixels = 512               # Pixels along x- and y- axis for each patch before extra pixels are added. x and y must have same number of pixels if processed with aug_hovernet.
    nr_pixels_overlap = 100     # Number of pixels overlap between patches
    # number of pixels along x- and y axis for each patch
    x = pixels + nr_extra_pixels
    y = pixels + nr_extra_pixels

    extractor = SlidingWindowPatchExtractor(
        img,
        patch_size=(x, y),
        stride=(x-nr_pixels_overlap, y-nr_pixels_overlap),
        pad_mode="constant",
        pad_constant_values=255,   # white padding
        within_bound=False
    )
    print(f"Number of images: {len(extractor)}")

    # Iterate over the patches and save each patch
    for patch_idx in range(len(extractor)): 
        patch = extractor[patch_idx]

        patch_image = Image.fromarray(patch)
        patch_filename = os.path.join(patch_save_path, f"patch_{patch_idx}.png")
        patch_image.save(patch_filename)


def _to_uint8_rgb(img: np.ndarray) -> np.ndarray:
    """Convert grayscale/RGB/RGBA image to uint8 RGB for visualization."""
    arr = np.asarray(img)

    if arr.ndim == 2:
        arr = np.stack([arr] * 3, axis=-1)
    elif arr.ndim == 3 and arr.shape[-1] == 4:
        arr = arr[..., :3]
    elif arr.ndim != 3 or arr.shape[-1] not in (3,):
        raise ValueError(f"Unsupported image shape for overview: {arr.shape}")

    if arr.dtype == np.uint8:
        return arr

    arr = arr.astype(np.float32)
    arr_min = arr.min()
    arr_max = arr.max()

    if arr_max == arr_min:
        return np.zeros_like(arr, dtype=np.uint8)

    arr = (arr - arr_min) / (arr_max - arr_min)
    arr = (arr * 255).clip(0, 255).astype(np.uint8)
    return arr


def extract_patches_and_create_overview_image(image_path: str | Path, patch_save_path: Path):
    patch_save_path = Path(patch_save_path)
    patch_save_path.mkdir(parents=True, exist_ok=True)

    try:
        img = tifffile.imread(image_path)
        image_name = Path(image_path).stem
    except Exception:
        img = image_path
        image_name = "image"

    nr_extra_pixels = 0       # Extra pixels per patch
    pixels = 4856              # Patch width/height before extra pixels
    nr_pixels_overlap = 0     # Overlap between neighboring patches

    patch_w = pixels + nr_extra_pixels
    patch_h = pixels + nr_extra_pixels
    stride_x = patch_w - nr_pixels_overlap
    stride_y = patch_h - nr_pixels_overlap

    extractor = SlidingWindowPatchExtractor(
        img,
        patch_size=(patch_w, patch_h),
        stride=(stride_x, stride_y),
        pad_mode="constant",
        pad_constant_values=255,   # white padding
        within_bound=False
    )

    print(f"Number of patches: {len(extractor)}")

    # ------------------------------------------------------------------
    # Save patches
    # ------------------------------------------------------------------
    for patch_idx in range(len(extractor)):
        patch = extractor[patch_idx]
        patch_image = Image.fromarray(patch)
        patch_filename = patch_save_path / f"patch_{patch_idx}.png"
        patch_image.save(patch_filename)

    # ------------------------------------------------------------------
    # Create overview image with patch numbers
    # ------------------------------------------------------------------
    h, w = img.shape[:2]

    # Assumes row-major patch ordering:
    # left-to-right across a row, then next row top-to-bottom.
    x_starts = list(range(0, w, stride_x)) or [0]
    y_starts = list(range(0, h, stride_y)) or [0]
    patch_positions = [(x0, y0) for y0 in y_starts for x0 in x_starts]

    if len(patch_positions) != len(extractor):
        print(
            f"Warning: computed {len(patch_positions)} positions but extractor has "
            f"{len(extractor)} patches. The overlay may be misaligned if the extractor "
            f"uses a different traversal order."
        )

    # Make a display-friendly RGB image
    overview_rgb = _to_uint8_rgb(img)

    # Downscale overview for easier viewing if needed
    overview_max_dim = 4000
    scale = min(1.0, overview_max_dim / max(h, w))
    overview_size = (int(round(w * scale)), int(round(h * scale)))

    overview_pil = Image.fromarray(overview_rgb).resize(
        overview_size, Image.Resampling.BILINEAR
    )
    draw = ImageDraw.Draw(overview_pil)

    # Font size based on patch size after scaling
    scaled_patch_w = max(1, int(round(patch_w * scale)))
    font_size = max(12, scaled_patch_w // 6)

    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", font_size)
    except Exception:
        font = ImageFont.load_default()

    for patch_idx, (x0, y0) in enumerate(patch_positions[:len(extractor)]):
        x1 = min(x0 + patch_w, w)
        y1 = min(y0 + patch_h, h)

        # Scale coordinates to overview image
        sx0 = int(round(x0 * scale))
        sy0 = int(round(y0 * scale))
        sx1 = int(round(x1 * scale))
        sy1 = int(round(y1 * scale))

        # Draw patch rectangle
        pistachio = (147, 197, 114)
        draw.rectangle([(sx0, sy0), (sx1, sy1)], outline=pistachio, width=2)

        # Draw patch index in center
        cx = (sx0 + sx1) // 2
        cy = (sy0 + sy1) // 2
        text = str(patch_idx)

        draw.text(
            (cx, cy),
            text,
            fill=pistachio,
            font=font,
            anchor="mm",
            stroke_width=2,
            stroke_fill=(255, 255, 255),
        )

    overview_path = patch_save_path / f"overview_image/{image_name}_patch_overview.png"
    if not overview_path.exists:
        overview_path.mkdir(parents=True)
    overview_pil.save(overview_path)
    print(f"Saved overview image to: {overview_path}")

def extract_patches_all_in_same_folder(root_dir: str | Path):
    root_dir = Path(root_dir)

    output_root = Path(
        "/media/jenny/Expansion/funcprost_Visium/20x/patches/second_group/"
    )

    # Add or remove extensions as needed
    image_extensions = {".tif"}

    image_paths = sorted(
        path
        for path in root_dir.iterdir()
        if path.is_file() and path.suffix.lower() in image_extensions
    )

    for image_path in image_paths:
        base_name = image_path.stem

        patch_save_path = output_root / base_name / "2048x2048"

        if patch_save_path.exists():
            print(
                f"Skipping {base_name}, output already exists: "
                f"{patch_save_path}"
            )
            continue

        patch_save_path.mkdir(parents=True, exist_ok=True)
        print(f"Directory {patch_save_path} was created")
        print(f"Processing: {image_path}")

        try:
            if image_path.suffix.lower() in {".tif", ".tiff"}:
                img = tifffile.imread(image_path)
            else:
                img = Image.open(image_path).convert("RGB")
                img = np.asarray(img)
        except Exception as error:
            print(f"Could not read {image_path}: {error}")
            continue

        pixels = 2048
        nr_extra_pixels = 0
        nr_pixels_overlap = 100

        x = pixels + nr_extra_pixels
        y = pixels + nr_extra_pixels

        extractor = SlidingWindowPatchExtractor(
            img,
            patch_size=(x, y),
            stride=(
                x - nr_pixels_overlap,
                y - nr_pixels_overlap,
            ),
            pad_mode="constant",
            pad_constant_values=255,
            within_bound=False,
        )

        print(f"Number of patches: {len(extractor)}")

        for patch_idx, patch in enumerate(extractor):
            patch_image = Image.fromarray(patch)
            patch_filename = patch_save_path / f"patch_{patch_idx}.png"
            patch_image.save(patch_filename)

        print(f"Finished processing {base_name}")

if __name__ == "__main__":

    # image_path = "/media/jenny/Expansion/MM_HE_masks/HE_MM089_C_270125_20x_BF_01/HE_MM089_C_270125_20x_BF_01_no_folds_mask.png"
    # image_path = "/media/.../patch_x.png"
    # image_path = "/media/jenny/Expansion/Prostata_Vilde/Analysis/Jenny/CD31/10x/images/Func015_CD31.tif"
    # patch_save_path = "/media/jenny/Expansion/Prostata_Vilde/Analysis/Jenny/CD31/10x/patches/patches_Func015/patches_Func015_16_squares/"
    
    # patch_save_path = Path("/media/jenny/Expansion/MM_HE_patches/HE_MM179_B_70225_20x_BF_01/aughovernet/2048x2048/")
    image_path = "/media/jenny/Expansion/HE_patches/20x/Func116_ST_HE_20x_BF_01/aughovernet/2048x2048/patch_24.png"
    # patch_mask_save_path = Path("/media/jenny/Expansion/MM_HE_patches/HE_MM172_2B_290125_20x_BF_01/aughovernet/2048x2048_mask/")
    patch_save_path = Path("/media/jenny/Expansion/test_nuclei/test3/patches_512x512/")
   
    # if not patch_mask_save_path.exists():
    #     patch_mask_save_path.mkdir(parents=True)
    #     print(f"Directory {patch_mask_save_path} was created")
    if not patch_save_path.exists():
        patch_save_path.mkdir(parents=True)
        print(f"Directory {patch_save_path} was created")
    
    # mask_path = "/media/.../masks/.../image.jpg"

    # extract_patches_wsi(wsi_path, patch_save_path)
    # extract_patches(image_path, patch_mask_save_path)
    extract_patches(image_path, patch_save_path)
    
    # extract_patches_all_in_folder(root_dir= Path("/media/jenny/Expansion/MetoxyLacc_HE_20x_TIFF/"))
    # extract_patches_all_in_folder_mask(root_dir= Path("/media/jenny/Expansion/MetoxyLacc_HE_20x_masks/"))

    # extract_patches_and_create_overview_image(image_path, patch_save_path)
    # extract_patches_all_in_same_folder(root_dir= "/media/jenny/Expansion/funcprost_Visium/20x/images/second_group/")