#!/usr/bin/env bash
# Builds the whole website into site/_site. Used by .github/workflows/pages.yml;
# can also be run locally from the repository root (needs Python and Ruby/Bundler).
set -euo pipefail
cd "$(dirname "$0")/.."

# 0. Stop with a clear message if something in content/ is wrong
python scripts/check_content.py

# File locations set in content/ (paths from the top folder of the repository)
cv_pdf="$(python scripts/content.py cv pdf)"

# 1. Gallery thumbnails and photo list, and the web version of the About-page photo
python scripts/build_gallery.py
python scripts/build_portrait.py

# 2. Date shown as "Last updated" on the Resume page: the last commit of the CV PDF
mkdir -p site/_data/generated
echo "cv_updated: $(git log -1 --format=%cs -- "$cv_pdf")" > site/_data/generated/build.yml

# 3. Jekyll pages, plus the public files from the repository root at the same
#    paths as in the repository (photos, CV PDF, profile.jpg, code folders, …;
#    see site/_plugins/root_files.rb).
#    SITE_URL / BASEURL (set by the workflow) publish the site somewhere else than
#    hezhouyi.github.io, e.g. a preview in a fork at yuying-jin.github.io/<repository>.
config=_config.yml
SITE_URL="$(printf '%s' "${SITE_URL:-}" | tr 'A-Z' 'a-z')"  # GitHub user names may have capitals; addresses are lower case
if [ -n "${SITE_URL:-}${BASEURL:-}" ]; then
  printf 'url: "%s"
baseurl: "%s"
' "${SITE_URL:-https://hezhouyi.github.io}" "${BASEURL:-}" > site/.build-override.yml
  config=_config.yml,.build-override.yml
fi
(cd site && bundle exec jekyll build --config "$config")

# 4. Title and author in the published CV's PDF metadata (for search results)
python scripts/set_pdf_metadata.py

echo "Site built in site/_site ($(du -sh site/_site | cut -f1))"
