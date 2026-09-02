import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv2u
import SimpleITK as sitk
from skimage.io import imread, imsave
from skimage.transform import rescale, rotate, AffineTransform, warp
from scipy.io import savemat

# Get the path to the project root directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# Append the project root to sys.path
if project_root not in sys.path:
    sys.path.append(project_root)
from tiatoolbox.wsicore.wsireader import WSIReader





# --------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------
def compute_scale_factor(fixed_image_rgb: np.ndarray, moving_image_rgb: np.ndarray) -> float:
    """
    Python equivalent of computeScaleFactor(fixed_image_rgb, moving_image_rgb) in MATLAB.
    """
    
    print("Read dimensions")
    fixed_h, fixed_w = fixed_image_rgb.shape[:2]
    moving_h, moving_w = moving_image_rgb.shape[:2]
    print("Calculate scale")
    scale_h = fixed_h / moving_h
    scale_w = fixed_w / moving_w
    print("Scale calculated")
    
    x = float((scale_h + scale_w) / 2.0)
    print("x is calculated")
    return x


def func_manual_rotation(moving_image_rgb: np.ndarray,
                         fixed_image_rgb: np.ndarray,
                         theta: float,
                         dx: float,
                         dy: float) -> np.ndarray:
    """
    Rough Python equivalent of funcManualRotation:
    rotate by theta (counterclockwise), then translate by (dx, dy).
    """
    rotated = rotate(
        moving_image_rgb,
        theta,
        resize=False,
        preserve_range=True,
    )

    # Translate: dx right, dy down
    tform = AffineTransform(translation=(dx, dy))
    warped = warp(rotated, tform.inverse, preserve_range=True)

    return warped.astype(moving_image_rgb.dtype)


def func_coreg(fixed_image_rgb: np.ndarray,
               moving_image_rgb: np.ndarray):
    """
    Python translation of MATLAB funcCoReg.

    MATLAB version (summary):
      - FileReduction = 0.2
      - fixed  = double(imresize(fixed_image_rgb(:,:,3), FileReduction));
      - moving = double(imresize(moving_image_rgb(:,:,1), FileReduction));
      - fixed  = abs(256-fixed); moving = abs(256-moving);
      - [optimizer,metric] = imregconfig('multimodal');
      - similarity, then affine registration
      - adjust translation by 1/FileReduction
      - imwarp each channel of moving_image_rgb with final tform
      - RotMatrix = tform.T

    Here we mimic that with SimpleITK. We downsample,
    invert intensities, register in 2 stages (similarity + affine),
    and apply the final transform to each channel.
    """
    FileReduction = 0.2

    # --- 1. Downsample and pick channels (3rd of fixed, 1st of moving) ---
    fixed_small = rescale(
        fixed_image_rgb[..., 2],     # fixed_image_rgb(:,:,3) in MATLAB (blue channel, 1-based)
        scale=FileReduction,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32)

    moving_small = rescale(
        moving_image_rgb[..., 0],    # moving_image_rgb(:,:,1) in MATLAB (red channel, 1-based)
        scale=FileReduction,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32)

    # --- 2. Invert intensities: abs(256 - image) ---
    fixed_small = np.abs(256.0 - fixed_small)
    moving_small = np.abs(256.0 - moving_small)

    # --- 3. Convert to SimpleITK images ---
    fixed_sitk = sitk.GetImageFromArray(fixed_small)
    moving_sitk = sitk.GetImageFromArray(moving_small)

    # Compensate for downsampling: keep same physical size as full-res images
    spacing_small = (1.0 / FileReduction, 1.0 / FileReduction)
    fixed_sitk.SetSpacing(spacing_small)
    moving_sitk.SetSpacing(spacing_small)

    # --- 4. Registration setup (multimodal-like) ---
    def configure_and_run_registration(initial_transform, fixed_img, moving_img):
        method = sitk.ImageRegistrationMethod()
        # Metric ~ imregconfig('multimodal')
        method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50) # “mutual-information metric” = a similarity measure that tells the optimizer how well two images align, by looking at how statistically dependent their intensities are at corresponding positions.
        method.SetMetricSamplingStrategy(method.RANDOM)
        method.SetMetricSamplingPercentage(0.1)

        method.SetInterpolator(sitk.sitkLinear)

        # Optimizer: approximate InitialRadius / MaximumIterations tweaks
        method.SetOptimizerAsRegularStepGradientDescent(
            learningRate=2.0,      # tune as needed
            minStep=1e-4,
            numberOfIterations=1000,
            relaxationFactor=0.5,
        )
        method.SetOptimizerScalesFromPhysicalShift()

        # Multi-resolution
        method.SetShrinkFactorsPerLevel(shrinkFactors=[4, 2, 1])
        method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2, 1, 0])
        method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

        method.SetInitialTransform(initial_transform, inPlace=False)

        final_transform = method.Execute(fixed_img, moving_img)
        return final_transform

    # --- 5. Stage 1: similarity transform (like imregtform(..., 'similarity', ...) ---
    #NOTE: Matlab performs no initial transformation unless it is specified by the user
    initial_similarity = sitk.CenteredTransformInitializer(
        fixed_sitk,
        moving_sitk,
        sitk.Similarity2DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY,
    )

    final_transform = configure_and_run_registration(
        initial_similarity,
        fixed_sitk,
        moving_sitk,
    )

    # --- 6. Stage 2: affine transform with similarity as initial guess ---

    # Needs to be converted from composite for the program to work
    if isinstance(similarity_transform, sitk.CompositeTransform):
        similarity_transform = similarity_transform.GetBackTransform()
        print("get back")
    else:
        similarity_transform = similarity_transform

    initial_affine = sitk.AffineTransform(2)
    initial_affine.SetCenter(similarity_transform.GetCenter())
    initial_affine.SetMatrix(similarity_transform.GetMatrix())
    initial_affine.SetTranslation(similarity_transform.GetTranslation())

    final_transform = configure_and_run_registration(
        initial_affine,
        fixed_sitk,
        moving_sitk,
    )

    # Needs to be converted from composite for the program to work
    if isinstance(final_transform, sitk.CompositeTransform):
        final_transform = final_transform.GetBackTransform()
    else: 
        final_transform = final_transform
    
 
    # NOTE:
    # In MATLAB you manually rescale tform.T(3,1:2) by 1/FileReduction.
    # Here we instead encoded the downsampling in the spacing, so the
    # transform is already defined in the "full-res" physical coordinate
    # system. No extra scaling is needed.

    # --- 7. Apply final transform to each channel of full-res moving image ---
    # Reference image: same size as fixed_image_rgb(:,:,3)
    ref_channel = fixed_image_rgb[..., 2].astype(np.float32)
    ref_sitk = sitk.GetImageFromArray(ref_channel)

    reg_channels = []
    for ch in range(moving_image_rgb.shape[2]):
        moving_ch = moving_image_rgb[..., ch].astype(np.float32)
        moving_ch_sitk = sitk.GetImageFromArray(moving_ch)

        reg_ch_sitk = sitk.Resample(
            moving_ch_sitk,
            ref_sitk,
            final_transform,
            sitk.sitkLinear,
            0.0,  # default fill
            moving_ch_sitk.GetPixelID(),
        )

        reg_ch = sitk.GetArrayFromImage(reg_ch_sitk)
        reg_channels.append(reg_ch)

    RegImage = np.stack(reg_channels, axis=-1).astype(moving_image_rgb.dtype)

    # --- 8. Build a 3x3 matrix similar to MATLAB's tform.T ---
    M = np.array(final_transform.GetMatrix()).reshape(2, 2)
    t = np.array(final_transform.GetTranslation())

    RotMatrix = np.array(
        [
            [M[0, 0], M[0, 1], 0.0],
            [M[1, 0], M[1, 1], 0.0],
            [t[0],    t[1],    1.0],
        ],
        dtype=float,
    )

    return RegImage, RotMatrix


# --------------------------------------------------------------------
# Main script logic (updated to use func_coreg)
# --------------------------------------------------------------------
name_list = ['179_A']

# Example choice of marker / stain
String = 'CD8'  # 'CD31', 'Ki67', 'PIN trippel', etc.

# Base directories
DirH = '/media/jenny/Expansion/MM_HE_pyramidal_tiff/Pyramidal_HE_MM179_A_70225_20x_BF_01.tif'
DirI = '/media/jenny/Expansion/co_reg/HE_MM179_A_70225_20x_BF_01/CD8/Pyramidal_Image_MM179_A_CD8.tif'  

# DirR = '/path/to/results/'         # where .mat files should be stored

# DirDigPat = os.path.join(
#     DirH,
#     String,
# )
# or:
# DirDigPat = os.path.join(DirH, 'Digital pathology', 'TIF-Files', String, '')

# Only process 3rd element (index 2), like: for i = 3 % numel(nameList)
for i in [0]:
    name = name_list[i]

    # fixed_path = os.path.join(DirI, f'{name} HE Visium.tif')
    # moving_path = os.path.join(DirDigPat, f'{name} {String}.tif')
    # Image paths
    moving_img_file_name = "/media/jenny/Expansion/MM_HE_pyramidal_tiff/Pyramidal_HE_MM009_B_270125.tif"
    fixed_img_file_name = "/media/jenny/Expansion/co_reg/HE_MM009_2_270125_20x_BF_01/CD8/Pyramidal_Image_MM009_B_CD8.tif"


    # Read images and reduce the resolution
    fixed_wsi_reader = WSIReader.open(input_img=fixed_img_file_name)
    # fixed_image_rgb = fixed_wsi_reader.slide_thumbnail(resolution=0.1563, units="power")
    fixed_image_rgb = fixed_wsi_reader.slide_thumbnail(resolution=0.2, units="power")
    moving_wsi_reader = WSIReader.open(input_img=moving_img_file_name)
    # moving_image_rgb = moving_wsi_reader.slide_thumbnail(resolution=0.1563, units="power")
    moving_image_rgb = moving_wsi_reader.slide_thumbnail(resolution=0.2, units="power")

    # moving_image_rgb = cv2.resize(moving_image_rgb, (fixed_image_rgb.shape[1], fixed_image_rgb.shape[0]), interpolation=cv2.INTER_CUBIC)
    fixed_image_rgb = cv2u.resize(fixed_image_rgb, (moving_image_rgb.shape[1], moving_image_rgb.shape[0]), interpolation=cv2u.INTER_CUBIC) # Makes sure fixed and moving image have same dimensions


    # # --- 1: Match the size of the Visium image ---
    # scale_factor = compute_scale_factor(fixed_image_rgb, moving_image_rgb)
    # moving_image_rgb_resized = rescale(
    #     moving_image_rgb,
    #     scale=scale_factor,
    #     channel_axis=-1,
    #     preserve_range=True,
    #     anti_aliasing=True,
    # ).astype(moving_image_rgb.dtype)

    # print("Scale factor calculated")

    # --- 2: Manual rotation (still hard-coded [20,0,0], like your temporary code) ---
    ManualRotation = np.array([0, 0, 0])  # [counterclockwise, right, down]
    Img_Rotate = func_manual_rotation(
        moving_image_rgb,
        fixed_image_rgb,
        ManualRotation[0],
        ManualRotation[1],
        ManualRotation[2],
    )

    # Evaluate rotation (rough imshowpair equivalent)
    plt.figure()
    # if fixed_image_rgb.ndim == 3 and Img_Rotate.ndim == 3:
    #     overlay = np.zeros_like(fixed_image_rgb, dtype=np.float32)
    #     overlay[..., 1] = fixed_image_rgb[..., 0]       # green
    #     overlay[..., 0] = Img_Rotate[..., 0]      # red + blue (magenta)
    #     overlay[..., 2] = Img_Rotate[..., 0]
    #     overlay = np.clip(overlay, 0, 255).astype(np.uint8)
    #     plt.imshow(overlay)
    # else:
    plt.imshow(fixed_image_rgb, alpha=0.8)
    plt.imshow(Img_Rotate, alpha=0.6)
    plt.title(f'Overlay (manual rot): {name}')
    plt.axis('off')

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    if fixed_image_rgb.ndim == 3:
        im0 = np.abs(256 - fixed_image_rgb[..., 2])
    else:
        im0 = np.abs(256 - fixed_image_rgb)
    axes[0].imshow(im0.astype(np.uint8), cmap='gray')
    axes[0].set_title('Fixed (proc)')
    axes[0].axis('off')

    if Img_Rotate.ndim == 3:
        im1 = np.abs(256 - Img_Rotate[..., 0])
    else:
        im1 = np.abs(256 - Img_Rotate)
    axes[1].imshow(im1.astype(np.uint8), cmap='gray')
    axes[1].set_title('Rotated (proc)')
    axes[1].axis('off')

    plt.tight_layout()
    plt.show()

    # --- 3: Co-register using translated func_coreg (funcCoReg) ---
    RegImage, RotMatrix = func_coreg(fixed_image_rgb, Img_Rotate)
    
    # --- 4: Visualize final registration (alpha-blended overlay) ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Make sure we have uint8 for display
    fixed_vis = fixed_image_rgb.astype(np.uint8)
    reg_vis   = RegImage.astype(np.uint8)

    axes[0].imshow(fixed_vis)
    axes[1].imshow(reg_vis, alpha=0.9) 

    plt.axis("off")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 8))

    # Make sure we have uint8 for display
    fixed_vis = fixed_image_rgb.astype(np.uint8)
    reg_vis   = RegImage.astype(np.uint8)

    plt.imshow(fixed_vis)
    plt.imshow(reg_vis, alpha=0.4) 

    plt.title("Final registration")
    plt.axis("off")
    plt.show()
    # # Save registered image
    # out_tif_path = os.path.join(DirI, f'{name} {String}.tif')
    # imsave(out_tif_path, RegImage.astype(np.uint8))

    # # Save transform + metadata (MATLAB-style .mat)
    # out_mat_path = os.path.join(DirR, f'{name} {String}.mat')
    # savemat(
    #     out_mat_path,
    #     {
    #         'RotMatrix': RotMatrix,
    #         'ManualRotation': ManualRotation,
    #         'scaleFactor': scale_factor,
    #     },
    # )