from pathlib import Path
import argparse


def rename_spaces(root_directory: str, dry_run: bool = False):
    root = Path(root_directory).resolve()

    if not root.exists():
        print(f"Directory does not exist: {root}")
        return

    if not root.is_dir():
        print(f"Not a directory: {root}")
        return

    # Work from the deepest items upward so renaming directories
    # doesn't interfere with paths we still need to process.
    items = sorted(
        root.rglob("*"),
        key=lambda path: len(path.parts),
        reverse=True
    )

    for item in items:
        if " " not in item.name:
            continue

        new_name = item.name.replace(" ", "_")
        new_path = item.with_name(new_name)

        # Avoid overwriting an existing file/directory.
        if new_path.exists():
            print(f"SKIPPED: {item} -> {new_path} (destination already exists)")
            continue

        print(f"{item} -> {new_path}")

        if not dry_run:
            item.rename(new_path)


if __name__ == "__main__":
    
    # Directory to process
    directory = Path("/media/jenny/Expansion/Prostata_Vilde/PIN_trippel/JPEG/")

    # Show what would be renamed without actually changing anything
    dry_run = False

    rename_spaces(directory, dry_run)