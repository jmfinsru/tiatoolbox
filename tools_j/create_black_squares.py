from PIL import Image
from pathlib import Path
import os


def generate_squares(output_dir: Path | str):
    # Store the images
    os.makedirs(output_dir, exist_ok=True)

    # Image size
    width, height = 2048, 2048

    # Create a single black image
    black_img = Image.new("RGB", (width, height), (0, 0, 0))

    # Generate x black images
    x = 1190
    for i in range(x):
        filename = f"patch_{i}.png"
        black_img.save(os.path.join(output_dir, filename))

    print("Created black images in:", output_dir)

if __name__ == "__main__":
    output_dir = Path("/media/.../")
    generate_squares(output_dir)