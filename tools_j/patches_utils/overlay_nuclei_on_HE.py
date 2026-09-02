from pathlib import Path
from PIL import Image
import numpy as np
import cv2

def replace_black_pixels(base_path: Path, input_path: Path, output_path: Path):
    base_img = Image.open(base_path).convert("RGBA")
    input_img = Image.open(input_path).convert("RGBA")

    if base_img.size != input_img.size:
        raise ValueError(
            f"Size mismatch: {base_path.name} is {base_img.size}, "
            f"but {input_path.name} is {input_img.size}"
        )

    base_arr = np.array(base_img)
    input_arr = np.array(input_img)

    # A pixel is considered black if R, G, and B are all 0.
    # Alpha is ignored when checking for black.
    black_mask = (
        (base_arr[:, :, 0] == 0) &
        (base_arr[:, :, 1] == 0) &
        (base_arr[:, :, 2] == 0)
    )

    # Replace black pixels in base image with matching pixels from input image
    base_arr[black_mask] = input_arr[black_mask]

    output_img = Image.fromarray(base_arr, mode="RGBA")
    output_img.save(output_path)


def create_nuclei_overlay(nuclei_mask_folder, HE_patches_folder, output_folder):
    print("start")
    output_folder.mkdir(parents=True, exist_ok=True)
    print(output_folder)
    base_images = sorted(nuclei_mask_folder.glob("overlay_patch_*.png"))
    #print(base_images)
    for base_path in base_images:
        # Extract number from overlay_patch_997.png
        number = base_path.stem.replace("overlay_patch_", "")

        input_path = HE_patches_folder / f"patch_{number}.png"
        output_path = output_folder / f"output_patch_{number}.png"

        if not input_path.exists():
            print(f"Skipping {base_path.name}: missing {input_path.name}")
            continue

        try:
            replace_black_pixels(base_path, input_path, output_path)
            print(f"Saved {output_path.name}")
        except Exception as e:
            print(f"Error processing {base_path.name}: {e}")



def create_nuclei_overlay_and_mid_region_marker(nuclei_mask_folder, HE_patches_folder, output_folder):
    print("start")
    output_folder.mkdir(parents=True, exist_ok=True)
    print(output_folder)

    base_images = sorted(nuclei_mask_folder.glob("overlay_patch_*.png"))

    for base_path in base_images:
        # Extract number from overlay_patch_997.png
        number = base_path.stem.replace("overlay_patch_", "")

        input_path = HE_patches_folder / f"patch_{number}.png"
        output_path = output_folder / f"output_patch_{number}.png"

        if not input_path.exists():
            print(f"Skipping {base_path.name}: missing {input_path.name}")
            continue

        try:
            # Create nuclei overlay
            replace_black_pixels(base_path, input_path, output_path)

            # Load created overlay
            image = cv2.imread(str(output_path))

            if image is None:
                print(f"Could not load {output_path.name}")
                continue

            height, width = image.shape[:2]

            # Draw rectangle marking the counted/owned region:
            # 50 pixels removed from every side
            border = 50

            cv2.rectangle(
                image,
                (border, border),
                (width - border - 1, height - border - 1),
                (0, 128, 255),   # Red in BGR
                2              # Line thickness
            )

            # Save image with rectangle
            cv2.imwrite(str(output_path), image)

            print(f"Saved {output_path.name}")

        except Exception as e:
            print(f"Error processing {base_path.name}: {e}")


if __name__ == "__main__":
    sample = "044"
    nuclei_mask_folder = Path(f"/media/jenny/Expansion/jenny_funcprost/pannuke/results/zhang_original_weights/correct_nuclei/Func{sample}_ST_HE_40x_BF_01/overlay/")
    HE_patches_folder = Path(f"/media/jenny/Expansion/HE_patches/40x/Func{sample}_ST_HE_40x_BF_01/aughovernet/2048x2048/")
    output_folder = Path(f"/media/jenny/Expansion/jenny_funcprost/pannuke/results/zhang_original_weights/correct_nuclei/Func{sample}_ST_HE_40x_BF_01/post_process_overlay/")
    create_nuclei_overlay(nuclei_mask_folder, HE_patches_folder, output_folder)