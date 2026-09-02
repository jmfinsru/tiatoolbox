from PIL import Image
from pathlib import Path
import os


def generate_squares(output_dir: Path | str):
    # Store the images
    os.makedirs(output_dir, exist_ok=True)

    # Image size
    width, height = 512, 512

    # Create a single black image
    black_img = Image.new("RGB", (width, height), (0, 0, 0))

    # Generate x black images
    x = 25
    for i in range(x):
        filename = f"patch_{i}.png"
        black_img.save(os.path.join(output_dir, filename))

    print("Created black images in:", output_dir)

if __name__ == "__main__":
    output_dir = Path("/media/jenny/Expansion/test_nuclei/test3/patches_512x512_masks/")
    generate_squares(output_dir)