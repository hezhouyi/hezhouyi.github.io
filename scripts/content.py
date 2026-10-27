"""Read the settings at the top of a file in content/ (between the two --- lines).

File locations in those settings (photo, photo folder, CV PDF, …) are written
from the top folder of the repository, starting with "/" and spelled as on
GitHub, e.g. "/profile.jpg"; that is also where the file appears on
the website. scripts/check_content.py checks them first.

Used by the other build scripts:   from content import setting
and by the shell scripts:          python scripts/content.py cv pdf
"""
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import yaml

ROOT = Path(__file__).resolve().parent.parent


def settings(name):
    text = (ROOT / "content" / f"{name}.md").read_text(encoding="utf-8")
    m = re.match(r"---\s*\n(.*?)\n---[ \t]*(\n|$)", text, re.S)
    return (yaml.safe_load(m.group(1)) or {}) if m else {}


def setting(name, key, required=True):
    """A file location from content/<name>.md, as a path inside the repository
    ("/Zhouyi%20He.pdf" -> "Zhouyi He.pdf")."""
    value = unquote(str(settings(name).get(key) or "").strip()).strip("/")
    if not value and required:
        sys.exit(f'error: content/{name}.md has no "{key}:" line')
    return value


if __name__ == "__main__":
    print(setting(*sys.argv[1:3]))
