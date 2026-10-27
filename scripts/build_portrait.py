"""Make the web version of the About-page photo, the file set as "photo:" in
content/about.md. Without a photo there, the About page shows none.

Any size and format works (JPG, PNG, WEBP, HEIC); the result is
site/assets/img/portrait.jpg, small enough for the web (3x the size it is
shown at), and its pixel size goes to site/_data/generated/portrait.json.
The original file is never changed.

Run from the repository root:  python scripts/build_portrait.py
"""
import json
import sys
from pathlib import Path

from PIL import Image, ImageOps

from content import setting

try:  # HEIC/HEIF support (iPhone photos)
    from pillow_heif import register_heif_opener
    register_heif_opener()
except ImportError:
    pass

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "site" / "assets" / "img" / "portrait.jpg"
DATA = ROOT / "site" / "_data" / "generated" / "portrait.json"
BOX = (792, 1056)  # the photo is shown 264 x 352 px


def main():
    photo = setting("about", "photo", required=False)
    if not photo:
        OUT.unlink(missing_ok=True)
        DATA.unlink(missing_ok=True)
        print("portrait: none set in content/about.md")
        return
    src = ROOT / photo
    if not src.is_file():
        sys.exit(f"error: photo {photo} (content/about.md) not found")

    Image.MAX_IMAGE_PIXELS = None
    with Image.open(src) as im:
        im = ImageOps.exif_transpose(im).convert("RGB")
        im.thumbnail(BOX, Image.LANCZOS)
        OUT.parent.mkdir(parents=True, exist_ok=True)
        im.save(OUT, "JPEG", quality=84, optimize=True, progressive=True)
        width, height = im.size

    DATA.parent.mkdir(parents=True, exist_ok=True)
    DATA.write_text(json.dumps({"width": width, "height": height}))
    print(f"portrait: {src.name} -> {width}x{height}, {OUT.stat().st_size // 1024} KB")


if __name__ == "__main__":
    main()
