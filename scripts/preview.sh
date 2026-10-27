#!/usr/bin/env bash
# Local preview without installing Ruby: builds the site inside Docker, then
# serves it at http://localhost:4000 (Ctrl+C to stop). Needs Docker and Python.
set -euo pipefail
cd "$(dirname "$0")/.."

# Docker on Windows needs a Windows path; `pwd -W` gives one in Git Bash.
repo="$(pwd -W 2>/dev/null || pwd)"

MSYS_NO_PATHCONV=1 docker run --rm -v "$repo:/repo" -w /repo \
  -e BUNDLE_PATH=vendor/bundle ruby:3.3 bash -c '
    set -e
    apt-get update -qq >/dev/null && apt-get install -y -qq python3-venv >/dev/null
    python3 -m venv /tmp/venv && . /tmp/venv/bin/activate
    pip install -q -r scripts/requirements.txt
    git config --global --add safe.directory /repo
    (cd site && bundle install --quiet)
    bash scripts/build_site.sh
  '

echo "Serving at http://localhost:4000"
python -m http.server 4000 --bind 127.0.0.1 --directory site/_site
