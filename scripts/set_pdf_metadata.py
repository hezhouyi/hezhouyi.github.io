"""Give the published CV a title and author, so search results show
"Zhouyi He — Curriculum Vitae" instead of the file name.

Only the copy in the built site (site/_site) is changed, and only fields that
are empty: once the uploaded PDF sets them (in LaTeX:
\\hypersetup{pdftitle=..., pdfauthor=...}), this script leaves them alone.

Run from the repository root after the site is built:  python scripts/set_pdf_metadata.py
"""
from pathlib import Path

from pypdf import PdfReader, PdfWriter

from content import setting

PDF = Path(__file__).resolve().parent.parent / "site" / "_site" / setting("cv", "pdf")
WANTED = {
    "/Title": "Zhouyi He — Curriculum Vitae",
    "/Author": "Zhouyi He",
    "/Subject": "Simulation and theory of biomolecular condensates",
}


def main():
    reader = PdfReader(PDF)
    current = reader.metadata or {}
    missing = {k: v for k, v in WANTED.items() if not str(current.get(k, "")).strip()}
    if not missing:
        print("CV metadata already set")
        return
    writer = PdfWriter(clone_from=reader)
    writer.add_metadata(missing)
    # Write to a new file first: the built copy may be a hard link to the original.
    tmp = PDF.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        writer.write(f)
    tmp.replace(PDF)
    print("CV metadata set:", ", ".join(k.lstrip("/") for k in missing))


if __name__ == "__main__":
    main()
