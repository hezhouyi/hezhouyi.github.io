"""Check the files in content/ before the site is built, and explain every
problem in plain words: which file, which line, and how to fix it.

If anything is wrong the build stops, so the website is not updated and the
previous version stays online. In GitHub Actions the problems are also shown
at the top of the run's page (Actions tab → the failed run).

Run from the repository root:  python scripts/check_content.py
"""
import os
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import yaml

ROOT = Path(__file__).resolve().parent.parent
CONTENT = ROOT / "content"
FILES = ["about", "code", "cv", "photography", "publications", "seo", "site"]

QUOTES = 'Put the whole text in "double quotes".'


class Settings(dict):
    """A block of settings that remembers the line of each of its keys."""
    line = 1
    lines = {}

    def line_of(self, key):
        return self.lines.get(key, self.line)


class Loader(yaml.SafeLoader):
    pass


def _settings(loader, node):
    s = Settings(loader.construct_mapping(node, deep=True))
    s.line = node.start_mark.line + 2  # + the opening --- line, and lines count from 1
    s.lines = {k.value: k.start_mark.line + 2 for k, _ in node.value}
    return s


Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _settings)


class Problems(list):
    def add(self, name, line, message, fix):
        self.append((f"content/{name}.md", line, message, fix))


def explain_yaml(e, lines):
    """(line, message, fix) for a YAML error, in words. YAML's own messages
    rarely name the real cause, so the text of the line is looked at too."""
    problem, context = e.problem or "", e.context or ""
    mark = e.problem_mark
    if "end of stream" in problem and e.context_mark:
        mark = e.context_mark
    i = min(mark.line, len(lines) - 1) if mark else 0
    text = lines[i] if lines else ""
    above = lines[i - 1] if i > 0 else ""
    value = text.split(":", 1)[1].strip() if ":" in text else text.strip()
    line = i + 2  # + the opening --- line, and lines count from 1
    unclosed = ("A quote mark is opened but not closed.", 'Add the closing " at the end of the text.')

    def quotes(s):
        return 0 if s.lstrip().startswith("#") else s.count('"')

    if "'\\t'" in problem:
        return line, "This line contains a tab.", "Use spaces instead of tabs."
    if quotes(text) % 2:
        return (line, *unclosed)
    if quotes(above) % 2:
        return (line - 1, *unclosed)
    if quotes(text) > 2:
        return line, 'The text contains a " inside the "double quotes" around it.', \
            "Inside the text use 'single quotes' or “curly quotes” instead."
    if value.startswith(("[", "{")) and not re.fullmatch(r"\[[^\]]*\]", value):
        return line, "This text starts with [ or {.", QUOTES
    if re.match(r"[@`*&!%]", value):
        return line, "This text starts with a special character (@ ` * & ! %).", QUOTES
    if "mapping values are not allowed" in problem and ": " in value:
        return line, 'This text contains a colon ":".', QUOTES
    return line, "This line does not line up with the lines around it.", \
        "Check that it starts with the same number of spaces as the lines like it (and the other lines of the same entry)."


def read(name, problems):
    path = CONTENT / f"{name}.md"
    if not path.is_file():
        problems.add(name, 1, "This file is missing.", "Restore it from the repository history (it is needed to build the site).")
        return None
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---"):
        return Settings()
    m = re.match(r"---[ \t]*\r?\n(.*?\r?\n)?---[ \t]*(\r?\n|$)", text, re.S)
    if not m:
        problems.add(name, 1, "The settings at the top start with a --- line, but there is no closing --- line.",
                     "Add a line with only --- after the last setting.")
        return None
    try:
        data = yaml.load(m.group(1) or "", Loader=Loader)
    except yaml.YAMLError as e:
        problems.add(name, *explain_yaml(e, (m.group(1) or "").splitlines()))
        return None
    if data is None:
        return Settings()
    if not isinstance(data, Settings):
        problems.add(name, 2, "The settings at the top are not written as  name: value  lines.",
                     "Write each setting on its own line, like   name: \"Zhouyi He\"")
        return None
    return data


def exists_exactly(rel):
    """Like Path.exists, but with the exact upper/lower case, as on GitHub."""
    here = ROOT
    for part in Path(rel).parts:
        if not here.is_dir() or part not in os.listdir(here):
            return False
        here = here / part
    return True


WEB = re.compile(r"(https?://|mailto:)", re.I)


def is_page(rel):
    """A page of the site itself (site/resume.html is served as /resume.html)."""
    return "/" not in rel and not rel.startswith(("_", ".")) and rel.endswith((".html", ".xml", ".txt")) \
        and exists_exactly(f"site/{rel}")


def location_problem(value, label, kind, example, suffix=None):
    """What is wrong with a link or file location, as (message, fix), or None.

    The one rule for content/: a web address (https://…, mailto:…), or "/" +
    the path of a file or page in this repository, spelled as on GitHub. kind is "file" or "folder" when
    it must be in this repository, "link" when a web address is fine too."""
    in_repo = f'"/" + the path of the {kind} in this repository, for example "{example}"'
    if WEB.match(value):
        if kind == "link":
            return None
        return (f"{label} is a web address, but it has to be a {kind} in this repository.",
                f"Upload the {kind} to the repository and write {in_repo}.")
    if not value.startswith("/"):
        if value.lower().startswith("www."):
            fix = f'Write "https://{value}".'
        elif exists_exactly(unquote(value).strip("/")) or (kind == "link" and is_page(value)):
            fix = f'Write "/{value}".'
        elif kind == "link":
            fix = f"Write a web address (https://…), or {in_repo}."
        else:
            fix = f"Write {in_repo}."
        return f'{label} "{value}" has to start with "/" (a file in this repository) or https:// (a website).', fix
    rel = unquote(re.split(r"[?#]", value, maxsplit=1)[0]).strip("/")
    if rel == "" and kind == "link":
        return None  # the home page
    what = "file or page" if kind == "link" else kind
    if not (exists_exactly(rel) or (kind == "link" and is_page(rel))):
        return (f'There is no {what} "/{rel}" in this repository.',
                "Check the spelling (upper and lower case count), or upload it there.")
    if kind == "link" and (ROOT / rel).is_dir() and not (ROOT / rel / "index.html").is_file():
        return (f'"/{rel}" is a folder, and a folder cannot be linked to.',
                "Link to a file in it, or to the folder's page on GitHub: open the folder on github.com and copy the address.")
    if kind in ("file", "folder") and (kind == "folder") != (ROOT / rel).is_dir():
        other = "folder" if kind == "file" else "file"
        return f'"/{rel}" is a {other}, not a {kind}.', f"Write {in_repo}."
    if suffix and not rel.lower().endswith(suffix):
        return f'"/{rel}" is not a {suffix} file.', f'Write a {suffix} file, for example "{example}".'
    return None


def check_location(problems, name, settings, key, example, kind="link", required=True, suffix=None):
    value = settings.get(key)
    line = settings.line_of(key)
    if value is None or str(value).strip() == "":
        if required:
            problems.add(name, line, f'The "{key}:" line is missing or empty.', f'Add it, for example   {key}: "{example}"')
        return
    problem = location_problem(str(value).strip(), f'"{key}:"', kind, example, suffix)
    if problem:
        problems.add(name, line, *problem)


def check_text_links(problems, name):
    """[text](link) inside any text follows the same rule."""
    path = CONTENT / f"{name}.md"
    if not path.is_file():
        return
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.lstrip().startswith("#"):  # comments (and headings)
            continue
        for target in re.findall(r"\]\(\s*<?([^)>]*?)>?\s*\)", line):
            target = re.sub(r'\s+"[^"]*"$', "", target.strip())  # [text](link "title")
            if target.startswith("#"):
                continue
            problem = location_problem(target, "The link", "link", "/resume.html")
            if problem:
                problems.add(name, n, *problem)


def entries(problems, name, settings, key, example):
    items = settings.get(key)
    if items is None:
        return []
    if not isinstance(items, list):
        problems.add(name, settings.line_of(key), f'The entries under "{key}:" must each start with "  - ".',
                     f"Compare with the example in the comment at the top of the file: {example}")
        return []
    good = []
    for item in items:
        if isinstance(item, Settings):
            good.append(item)
        else:
            problems.add(name, settings.line_of(key), f'One of the entries under "{key}:" is not written as  name: value  lines.',
                         f"Compare with the other entries: {example}")
    return good


def require(problems, name, item, keys, what):
    for key in keys:
        if item.get(key) in (None, ""):
            problems.add(name, item.line, f'This {what} has no "{key}:".', f'Add a line   {key}: "…"   to it.')


def main():
    sys.stdout.reconfigure(encoding="utf-8")  # for → and quotes in a Windows console
    problems = Problems()
    data = {name: read(name, problems) for name in FILES}

    for name in FILES:
        if data[name] is not None:
            check_text_links(problems, name)

    if (about := data["about"]) is not None:
        check_location(problems, "about", about, "photo", "/profile.jpg", kind="file", required=False)
        for b in entries(problems, "about", about, "buttons", 'label: "ORCID"  link: "https://orcid.org/…"'):
            check_location(problems, "about", b, "link", "/resume.html")
        for e in entries(problems, "about", about, "elsewhere", 'label: "Blog"  qrcode: "/content/qrcodes/qr-wechat.png"'):
            check_location(problems, "about", e, "qrcode", "/content/qrcodes/qr-wechat.png", kind="file", required=False)
            check_location(problems, "about", e, "link", "https://space.bilibili.com/…", required=False)

    if (site := data["site"]) is not None:
        for l in entries(problems, "site", site, "footer_links", 'label: "ORCID"  link: "https://orcid.org/…"'):
            check_location(problems, "site", l, "link", "https://orcid.org/…")

    if (seo := data["seo"]) is not None and isinstance(person := seo.get("person"), Settings):
        for key in ("works_at", "studied_at"):
            for o in entries(problems, "seo", person, key, 'name: "TU Dresden"  url: "https://tu-dresden.de/"'):
                check_location(problems, "seo", o, "url", "https://tu-dresden.de/", required=False)
        profiles = person.get("profiles") or []
        for url in profiles if isinstance(profiles, list) else [profiles]:
            if not re.match(r"https?://", str(url)):
                problems.add("seo", person.line_of("profiles"), f'The profile "{url}" is not a web address.',
                             'Write the full address, starting with https://')

    if (cv := data["cv"]) is not None:
        check_location(problems, "cv", cv, "pdf", "/Zhouyi He.pdf", kind="file", suffix=".pdf")

    if (photography := data["photography"]) is not None:
        check_location(problems, "photography", photography, "folder", "/image", kind="folder")

    code_ids = set()
    if (code := data["code"]) is not None:
        for c in entries(problems, "code", code, "code", "id: polymer-lattice  name: \"…\"  url: \"https://…\""):
            require(problems, "code", c, ["id", "name"], "code entry")
            check_location(problems, "code", c, "url", "https://github.com/hezhouyi/hezhouyi.github.io/tree/master/…", required=False)
            check_location(problems, "code", c, "docs", "/Polymer lattice simulation/manual.pdf", required=False)
            if c.get("id") in code_ids:
                problems.add("code", c.line, f'Two code entries have the id "{c["id"]}".', "Give each one its own id.")
            code_ids.add(c.get("id"))

    if (pubs := data["publications"]) is not None:
        for p in entries(problems, "publications", pubs, "publications", 'title: "…"  authors: "…"  venue: "…"  year: 2026'):
            require(problems, "publications", p, ["title", "authors", "venue", "year"], "paper")
            check_location(problems, "publications", p, "pdf", "/papers/He 2025.pdf", required=False)
            check_location(problems, "publications", p, "data", "https://doi.org/…", required=False)
            if p.get("year") not in (None, "") and not re.fullmatch(r"\d{4}", str(p["year"])):
                problems.add("publications", p.line_of("year"), f'"year: {p["year"]}" is not a year.', "Write only the year, for example   year: 2026")
            ids = p.get("code")
            if ids is None:
                continue
            if not isinstance(ids, list):
                problems.add("publications", p.line_of("code"), '"code:" must be a list in square brackets.',
                             f"Write it like   code: [{ids}]")
                continue
            if data["code"] is None:
                continue
            for i in ids:
                if i not in code_ids:
                    known = ", ".join(sorted(map(str, code_ids))) or "none yet"
                    problems.add("publications", p.line_of("code"), f'"code:" mentions "{i}", which is not an id in content/code.md.',
                                 f"Use one of the ids there ({known}), or add a code entry with this id.")

    if not problems:
        print("content/: everything looks fine")
        return 0
    report(problems)
    return 1


def report(problems):
    n = len(problems)
    head = (f"{n} problem{'s' if n > 1 else ''} in content/ — the website was NOT updated; "
            "the previous version stays online until this is fixed.")
    print(head + "\n")
    for file, line, message, fix in problems:
        print(f"  {file}, line {line}\n    {message}\n    → {fix}\n")

    if os.environ.get("GITHUB_ACTIONS"):
        for file, line, message, fix in problems:
            text = f"{message} → {fix}".replace("%", "%25").replace("\n", "%0A")
            print(f"::error file={file},line={line},title={file} line {line}::{text}")
        summary = os.environ.get("GITHUB_STEP_SUMMARY")
        if summary:
            with open(summary, "a", encoding="utf-8") as f:
                f.write(f"### {head}\n\n")
                for file, line, message, fix in problems:
                    f.write(f"- **{file}, line {line}**: {message}<br>→ {fix}\n")
                f.write("\nEdit the file on GitHub (pencil icon), fix the line, and commit; the site then rebuilds by itself.\n")


if __name__ == "__main__":
    sys.exit(main())
