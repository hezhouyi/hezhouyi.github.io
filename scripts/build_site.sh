#!/usr/bin/env bash
# Builds the whole website into site/_site. Used by .github/workflows/pages.yml;
# can also be run locally from the repository root (needs Python and Ruby/Bundler).
set -euo pipefail
cd "$(dirname "$0")/.."

# 0. Stop with a clear message if something in content/ is wrong
python scripts/check_content.py

# 1. Gallery thumbnails and photo list, and the web version of the About-page photo
python scripts/build_gallery.py
python scripts/build_portrait.py

# 2. Jekyll pages, plus the public files from the repository root at the same
#    paths as in the repository (photos, CV PDF, profile.jpg, code folders, …;
#    see site/_plugins/root_files.rb). The Resume page's "Last updated" date is
#    looked up in git by site/_plugins/content.rb.
#    SITE_URL / BASEURL (set by the workflow) publish the site somewhere else than
#    hezhouyi.github.io, e.g. a preview in a fork at yuying-jin.github.io/<repository>.
config=_config.yml
SITE_URL="$(printf '%s' "${SITE_URL:-}" | tr 'A-Z' 'a-z')"  # GitHub user names may have capitals; addresses are lower case
if [ -n "${SITE_URL:-}${BASEURL:-}" ]; then
  printf 'url: "%s"\nbaseurl: "%s"\n' "${SITE_URL:-https://hezhouyi.github.io}" "${BASEURL:-}" > site/.build-override.yml
  config=_config.yml,.build-override.yml
fi
(cd site && bundle exec jekyll build --config "$config")

# 3. Title and author in the published CV's PDF metadata (for search results)
python scripts/set_pdf_metadata.py

echo "Site built in site/_site ($(du -sh site/_site | cut -f1))"
