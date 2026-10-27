"""Make the Photography page's images from the photos in the folder set in
content/photography.md ("folder:", normally image/).

For every photo there this writes
  site/assets/gallery/thumb/<name>.jpg   800 px on the long side, for the grid
  site/assets/gallery/large/<name>.jpg  2000 px on the long side, for the lightbox
and lists them, with caption and size, in site/_data/generated/photos.json.

The originals are never changed. Outputs are reused while the
original's size is unchanged (the deploy workflow caches assets/gallery/), and
removed when the original is deleted.

Run from the repository root:  python scripts/build_gallery.py
"""
import json
import re
import sys
from pathlib import Path
from urllib.parse import quote

from PIL import Image, ImageOps

from content import setting

try:  # HEIC/HEIF support (iPhone photos)
    from pillow_heif import register_heif_opener
    register_heif_opener()
    HEIF = True
except ImportError:
    HEIF = False

ROOT = Path(__file__).resolve().parent.parent
FOLDER = setting("photography", "folder").rstrip("/")
SOURCE = ROOT / FOLDER
OUT = ROOT / "site" / "assets" / "gallery"
DATA = ROOT / "site" / "_data" / "generated" / "photos.json"
MANIFEST = OUT / "manifest.json"

SIZES = {"thumb": (800, 78), "large": (2000, 85)}  # long side in px, JPEG quality
EXTENSIONS = {".jpg", ".jpeg", ".png", ".webp", ".gif", ".bmp", ".tif", ".tiff"}
if HEIF:
    EXTENSIONS |= {".heic", ".heif"}

Image.MAX_IMAGE_PIXELS = None  # some originals are very large camera files


def caption(filename):
    """'DresdenFrauenKirche.jpg' -> 'Dresden Frauen Kirche', 'Oslo1.jpg' -> 'Oslo'."""
    s = Path(filename).stem
    s = re.sub(r"[_.]+", " ", s)
    s = re.sub(r"([a-z])([A-Z])", r"\1 \2", s)
    s = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1 \2", s)
    s = re.sub(r"\s*\d+$", "", s)
    return re.sub(r"\s+", " ", s).strip()


def output_name(filename):
    # Keep the original extension in the name so "Oslo.jpg" and "Oslo.png" cannot collide.
    return filename + ".jpg"


def render(src, name):
    with Image.open(src) as im:
        # Decode JPEGs at reduced size: much faster for 50 MB camera files.
        im.draft("RGB", (SIZES["large"][0], SIZES["large"][0]))
        im = ImageOps.exif_transpose(im).convert("RGB")
        sizes = {}
        for kind, (side, quality) in SIZES.items():
            copy = im.copy()
            copy.thumbnail((side, side), Image.LANCZOS)
            copy.save(OUT / kind / name, "JPEG", quality=quality, optimize=True, progressive=True)
            sizes[kind] = copy.size
    return sizes["thumb"]


def main():
    if not SOURCE.is_dir():
        sys.exit(f"error: photo folder {FOLDER}/ (content/photography.md) not found")
    for kind in SIZES:
        (OUT / kind).mkdir(parents=True, exist_ok=True)
    DATA.parent.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}

    photos, made, reused, failed = [], 0, 0, []
    files = sorted((p for p in SOURCE.iterdir() if p.is_file() and p.suffix.lower() in EXTENSIONS),
                   key=lambda p: p.name.lower())
    for src in files:
        name = output_name(src.name)
        size = src.stat().st_size
        entry = manifest.get(src.name)
        outputs_exist = all((OUT / kind / name).exists() for kind in SIZES)
        try:
            if entry and entry["bytes"] == size and outputs_exist:
                width, height = entry["width"], entry["height"]
                reused += 1
            else:
                width, height = render(src, name)
                made += 1
        except Exception as e:  # one unreadable photo must not break the site
            failed.append(f"{src.name}: {e}")
            continue
        manifest[src.name] = {"bytes": size, "width": width, "height": height}
        photos.append({
            "file": src.name,
            "caption": caption(src.name),
            "original": "/" + quote(f"{FOLDER}/{src.name}"),
            "thumb": "/assets/gallery/thumb/" + quote(name),
            "large": "/assets/gallery/large/" + quote(name),
            "width": width,
            "height": height,
        })

    # Forget photos that were deleted from the folder
    current = {src.name for src in files}
    for old in [k for k in manifest if k not in current]:
        for kind in SIZES:
            (OUT / kind / output_name(old)).unlink(missing_ok=True)
        del manifest[old]

    MANIFEST.write_text(json.dumps(manifest, indent=1))
    DATA.write_text(json.dumps(photos, indent=1, ensure_ascii=False), encoding="utf-8")

    print(f"{len(photos)} photos: {made} processed, {reused} reused from cache")
    if not HEIF:
        for p in SOURCE.iterdir():
            if p.suffix.lower() in {".heic", ".heif"}:
                print(f"warning: skipped {p.name} (install pillow-heif for HEIC support)")
    for f in failed:
        print(f"warning: could not read {f}")


if __name__ == "__main__":
    sys.exit(main())
