import numpy as np
import os
from skimage import io, color, morphology, filters, measure, segmentation
from scipy import ndimage as ndi
from tifffile import imwrite


def create_tumor_mask(
    rgb,
    min_tissue_area=50_000,
    min_hole_area=5_000,
):
    """
    Create tumor/tissue ROI mask from RGB histology image.

    Parameters
    ----------
    rgb : ndarray (H, W, 3)
        RGB histology image
    min_tissue_area : int
        Minimum connected tissue area to keep
    min_hole_area : int
        Fill holes smaller than this size

    Returns
    -------
    tumor_mask : ndarray (H, W), bool
    """

    # Normalize image safely
    rgb = rgb.astype(np.float64)
    rgb /= rgb.max()

    # Convert to optical density
    od = -np.log(rgb + 1e-8)

    # Tissue vs background (OD-based)
    # Background ≈ very low OD across channels
    tissue_mask = np.mean(od, axis=2) > 0.05

    tissue_mask = morphology.binary_opening(
        tissue_mask, morphology.disk(5)
    )

    tissue_mask = morphology.binary_closing(
        tissue_mask, morphology.disk(10)
    )

    tissue_mask = ndi.binary_fill_holes(tissue_mask)

    # Hematoxylin enrichment (nuclei-dense regions)
    hed = color.separate_stains(rgb, color.hed_from_rgb)
    hematoxylin = hed[..., 0]

    # Adaptive threshold (robust across slides)
    h_thresh = filters.threshold_otsu(hematoxylin[tissue_mask])
    nuclei_rich = hematoxylin > h_thresh

    nuclei_rich &= tissue_mask

    nuclei_rich = morphology.binary_closing(
        nuclei_rich, morphology.disk(12)
    )

    nuclei_rich = ndi.binary_fill_holes(nuclei_rich)

    # Remove small regions
    labeled = measure.label(nuclei_rich)
    tumor_mask = np.zeros_like(nuclei_rich, dtype=bool)

    for region in measure.regionprops(labeled):
        if region.area >= min_tissue_area:
            tumor_mask[labeled == region.label] = True

    # Final cleanup
    tumor_mask = morphology.remove_small_holes(
        tumor_mask, min_hole_area
    )

    tumor_mask = morphology.binary_closing(
        tumor_mask, morphology.disk(20)
    )

    tumor_mask = ndi.binary_fill_holes(tumor_mask)

    return tumor_mask


# Helper functions
def area_filter(mask, min_area, max_area):
    labeled = measure.label(mask)
    props = measure.regionprops(labeled)
    keep = np.zeros_like(mask, dtype=bool)
    for p in props:
        if min_area <= p.area <= max_area:
            keep[labeled == p.label] = True
    return keep


def eccentricity_filter(mask, max_ecc=0.95):
    labeled = measure.label(mask)
    keep = np.zeros_like(mask, dtype=bool)
    for p in measure.regionprops(labeled):
        if p.eccentricity <= max_ecc:
            keep[labeled == p.label] = True
    return keep

def apply_tumor_mask(namelist, input_path_for_tissue_mask):
    
    for name in namelist:
        img = io.imread(os.path.join(input_path_for_tissue_mask, f'{name}.png'))
        print(f"Creating tumor mask for {name}")
        print("------------------------------------")
        tumor_mask = create_tumor_mask(img)
        imwrite(os.path.join(save_path, f'{name}_tumor_mask.png'),
                tumor_mask.astype(np.uint8) * 255)

def ki67(namelist, save_path, input_path_color_deconvolution):

    # Proliferating cells (Ki67)
    IL = 0.12
    se8 = morphology.disk(8)

    for name in namelist:
        brownImg = io.imread(os.path.join(input_path_color_deconvolution, f'{name}_dab_od.png')).astype(float)
        tumor_mask = io.imread(os.path.join(save_path, f'{name}_tumor_mask.png')) > 0
        print(f"Performing morphological operations for {name}")
        print("------------------------------------")
        prolif = brownImg > IL
         
        prolif = morphology.binary_dilation(prolif, se8)
        prolif = morphology.binary_erosion(prolif, se8)
        prolif = ndi.binary_fill_holes(prolif)

        prolif = area_filter(prolif, 70, 2000)
        prolif &= tumor_mask
        
        # Binary mask for proliferating cells
        imwrite(os.path.join(save_path, f'{name}_BW_Ki67.png'),
                prolif.astype(np.uint8) * 255)

def hematoxylin(namelist, input_path_for_tissue_mask, save_path):
    # Hematoxylin cell mask
    backThresh = 80
    se8 = morphology.disk(20)
    se3 = morphology.disk(8)

    for name in namelist:
        rgb = io.imread(os.path.join(input_path_for_tissue_mask, f'{name}.png'))
        BW_Ki67 = io.imread(os.path.join(save_path, f'{name}_BW_Ki67.png')) > 0
        tumor_mask = io.imread(os.path.join(save_path, f'{name}_tumor_mask.png')) > 0
        print("1")
        target = rgb[..., 0]
        print("2")
        prolif_inv = ~morphology.binary_dilation(BW_Ki67, se3)
        print("3")
        cellMask = target <= backThresh
        cellMask &= prolif_inv & tumor_mask
        print("4")
        cellMask = ndi.binary_fill_holes(cellMask)
        print("5")
        # -- Lower threshold for large objects
        large = area_filter(cellMask, 1000, np.inf)
        print("6")
        largeMask = (target <= 50) & large
        print("7")
        cellMask = (cellMask & ~large) | largeMask
        print("8")
        # -- Merge small fragments
        small = cellMask & ~morphology.remove_small_objects(cellMask, min_size=51)
        print("9")
        small = morphology.binary_dilation(small, se8)
        print("10")
        small = morphology.binary_erosion(small, se8)
        print("11")
        small = ndi.binary_fill_holes(small)
        print("12")

        cellMask |= small
        print(f"Performing watershed for {name}")
        print("------------------------------------")
        # -- Watershed for large nuclei
        large = area_filter(cellMask, 350, np.inf)
        print("13")
        D = ndi.distance_transform_edt(large)
        print("14")
        labels = segmentation.watershed(-D, markers=None, mask=large)
        print(len(labels))
        cellMask = (cellMask & ~large) | (labels > 0)
        print("15")

        cellMask = ndi.binary_fill_holes(cellMask)
        print("16")
        # -- Final filters
        cellMask = area_filter(cellMask, 40, 700)
        print("17")
        cellMask = eccentricity_filter(cellMask, max_ecc=0.95)
        print(len(cellMask))
        imwrite(os.path.join(save_path, f'{name}_BW_CellMask_BackThres_ny2.png'),
                cellMask.astype(np.uint8) * 255)

        print('Pipeline complete')

if __name__ == "__main__":
    
    # Paths
    input_path_for_tissue_mask = '/media/jenny/Expansion/Almac_TMA_ki67/T51-001/T51-001_Ki67_1_XR_40-2025-11-25_11.06.47/full/'
    input_path_color_deconvolution = '/media/jenny/Expansion/Almac_TMA_ki67_results/T51-001/'
    save_path = '/media/jenny/Expansion/Almac_TMA_ki67_results/T51-001/'

    namelist = ['T51-001_Ki67_1_XR_40-2025-11-25_11.06.47_A0a']
    os.makedirs(save_path, exist_ok=True)
    
    # apply_tumor_mask()
    # ki67(namelist=namelist, save_path=save_path, input_path_color_deconvolution=input_path_color_deconvolution)
    hematoxylin(namelist=namelist, input_path_for_tissue_mask=input_path_for_tissue_mask, save_path=save_path)