import os
import numpy as np
from skimage import io, color, img_as_float, exposure
from tifffile import imwrite


# Paths
input_path = '/media/jenny/Expansion/Almac_TMA_ki67/T51-001/T51-001_Ki67_1_XR_40-2025-11-25_11.06.47/full/'
save_path = '/media/jenny/Expansion/Almac_TMA_ki67_results/T51-001/'

nameList = ['T51-001_Ki67_1_XR_40-2025-11-25_11.06.47_A0a']

# Choose stain matrix
# HED = Hematoxylin, Eosin, DAB
stain_matrix = color.hed_from_rgb
print(stain_matrix)
# Process samples
for name in nameList:
    print(f'Processing {name}')

    file_path = os.path.join(
        input_path,
        f'{name}.png'   
    )

    # Load image (float in [0,1])
    rgb = img_as_float(io.imread(file_path))
    
    # Color deconvolution
    
    # Output is optical density (OD) space
    # Shape: (H, W, 3)
    hed = color.separate_stains(rgb, stain_matrix)

    # Channels:
    # hed[..., 0] → Hematoxylin
    # hed[..., 1] → Eosin
    # hed[..., 2] → DAB (or background if no DAB)

    dab_od = hed[..., 2]
    
    # Save as float32 png
    imwrite(
        os.path.join(save_path, f'{name}_dab_od.png'),
        dab_od.astype(np.float32)
    )

    # Shift to make minimum zero
    dab_shifted = dab_od - np.min(dab_od)

    # Rescale to 0–1 for visualization
    dab_vis = exposure.rescale_intensity(dab_shifted, out_range=(0,1))

     # Save as float32 png
    imwrite(
        os.path.join(save_path, f'{name}_dab_vis.png'),
        dab_vis.astype(np.float32)
    )

print('Color deconvolution finished')

