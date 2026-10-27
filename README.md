# hezhouyi.github.io

Hi, this is Zhouyi He's personal github, it contains my personal webpage and research code and data associated with publications.

The site is built with [Jekyll](https://jekyllrb.com/) by a GitHub Actions workflow
(`.github/workflows/pages.yml`). Every change pushed to `master` is published
automatically within a few minutes — nothing needs to be run on your computer.

## Updating the site

All of these can be done in the GitHub web editor (open the file, click the pencil, commit).

| To change… | Edit |
|---|---|
| About page: title, card, buttons, greeting, About text, Interests, Elsewhere | `content/about.md` |
| Name at the top of every page (with 贺周易), links at the bottom | `content/site.md` |
| Titles and descriptions in search results, link previews, the profile for search engines | `content/seo.md` |
| About-page photo | Upload the photo anywhere in the repository and set `photo:` in `content/about.md` to it, now `/profile.jpg` (any size; resized automatically) |
| Publications | `content/publications.md` — newest first |
| Code shown on the Publications & Code page, and which paper uses it | `content/code.md` and the `code:` field of each paper |
| Photos | Upload to or delete from `image/` (the folder is set by `folder:` in `content/photography.md`). Thumbnails are made automatically; the caption comes from the file name (`DresdenFrauenKirche.jpg` → "Dresden Frauen Kirche") |
| CV | Upload the new PDF on GitHub ("Add file" → "Upload files") with the same name, now "Zhouyi He.pdf". Its location is set by `pdf:` in `content/cv.md`; "Last updated" follows the date of the upload |

**Links and file locations** in `content/` are always one of two kinds:

- a full web address: `https://…` (or `mailto:…`);
- a file, folder or page of this site, written from the top folder of the repository and
  starting with `/`, spelled exactly as on GitHub (spaces allowed, upper/lower case matters):
  `/profile.jpg`, `/Zhouyi He.pdf`, `/image`, `/resume.html`.

If something is wrong (a missing file, a typo in the settings), the build stops with a
message saying which file and line to fix, and the site online stays as it was. The
message is at the top of the failed run in the repository's **Actions** tab.

## How the site is built

Updating the site never requires installing anything: GitHub Actions does all of this.

- `content/` holds everything you edit: one Markdown file per page or list (settings at
  the top, text below), and the QR codes in `content/qrcodes/`. `scripts/check_content.py` checks them before every build.
- `site/` holds the Jekyll templates, layouts and styles that turn `content/` into pages.
- Public files at the repository root (`image/`, `Zhouyi He.pdf`, `profile.jpg`, code
  folders, …) are published next to the pages at the same paths, so a file's link
  changes only if the file is moved or renamed (`site/_plugins/root_files.rb`).
- `scripts/build_site.sh` runs the whole build: it makes the gallery thumbnails with
  Python (Jekyll cannot resize images), then runs Jekyll. For a quick preview after
  that, `cd site && bundle exec jekyll serve` shows the whole site, files included, and
  rebuilds it when anything in `content/` or a public file changes (photos and the
  About-page photo need `scripts/build_gallery.py` / `scripts/build_portrait.py` first).

To preview the full site on your own computer, with only Docker and Python
installed (no Ruby needed), run from the repository root:

```sh
bash scripts/preview.sh
```

It builds the site inside Docker and opens it at http://localhost:4000
(Ctrl+C to stop). Run it again after changing anything to see the update.
