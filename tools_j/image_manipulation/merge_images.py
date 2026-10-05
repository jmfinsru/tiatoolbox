from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


def combine_images(folder, base_name, I_ref, extension=".png"):
    folder = Path(folder)

    image_paths = [
        folder / f"{base_name}{i}{extension}"
        for i in I_ref
    ]

    # Check that all images exist
    for path in image_paths:
        if not path.exists():
            raise FileNotFoundError(f"Could not find: {path}")

    images = [Image.open(path).convert("RGB") for path in image_paths]
    
    font = font = ImageFont.truetype("DejaVuSans.ttf", 40)

    label_height = 60
    
    total_width = sum(img.width for img in images)
    max_height = max(img.height for img in images)

    combined = Image.new(
        "RGB",
        (total_width, max_height + label_height),
        "white"
    )

    draw = ImageDraw.Draw(combined)

    x = 0

    for img, ref_value in zip(images, I_ref):
        # Center image vertically within image area
        y = max_height - img.height

        combined.paste(img, (x, y))

        # Label text
        label = f"I_ref = {ref_value}"

        # Find text dimensions
        bbox = draw.textbbox((0, 0), label, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]

        # Center label under image
        text_x = x + (img.width - text_width) // 2
        text_y = max_height + (label_height - text_height) // 2

        draw.text(
            (text_x, text_y),
            label,
            fill="black",
            font=font
        )

        x += img.width

    output_path = folder / f"{base_name}combined.png"
    combined.save(output_path)

    print(f"Saved to: {output_path}")



if __name__ == "__main__":
    region_mode = "square" # circle or square, where square includes more tissue around the visium spots while the circle only look at the visium spots
    measurement_type = "pixels" # cells or pixels
    I_ref = ["0.5", "0.7", "1.0"]
    column_name = "severity_adjusted_AMACR_pixel"
    downsample_factor = 10
    biopsies = ["Func006","Func006_LN", "Func009", "Func015", "Func028", "Func029", "Func033", "Func036", "Func043", "Func044", "Func050", "Func050_LN", "Func052", "Func082", "Func083", "Func089", "Func097", "Func112", "Func116", "Func117"]
    # biopsies = ["Func117"]
    
    for biopsy in biopsies:

        input_folder = f"/media/jenny/Expansion/Prostata_Vilde/Co-reg_images_10x/results/results_amacr_{measurement_type}_in_spots/{biopsy}/{region_mode}/"
        base_name = f"{column_name}_overlay_heatmap_I_ref_"

        combine_images(
            folder=input_folder,
            base_name=base_name,
            I_ref=I_ref,
            extension=".png",
        )