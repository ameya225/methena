#!/usr/bin/env bash
# Publish the methena repository to GitHub.
# Usage (run inside the unpacked methena/ folder):
#   ./publish_to_github.sh <github-user-or-org> [repo-name] [--private]
# Requires git. If the GitHub CLI 'gh' is installed and logged in (gh auth login), the repository
# is created and pushed automatically; otherwise the script prints the two commands to run.
set -euo pipefail
OWNER="${1:?Usage: $0 <github-user-or-org> [repo-name] [--private]}"
NAME="${2:-methena}"
VIS="--public"; [[ "${3:-}" == "--private" || "${2:-}" == "--private" ]] && VIS="--private"
[[ "$NAME" == "--private" ]] && NAME=methena
[ -f DESCRIPTION ] && [ -d .git ] || { echo "Run this inside the methena repository folder."; exit 1; }
# fill in the owner/repo placeholders in README badge/links, DESCRIPTION and CITATION.cff
for f in README.md DESCRIPTION CITATION.cff; do
  sed -i.bak "s#OWNER/methena#${OWNER}/${NAME}#g" "$f" && rm -f "$f.bak"
done
git add -A
git -c user.name="${GIT_AUTHOR_NAME:-$(git config user.name || echo methena)}" \
    -c user.email="${GIT_AUTHOR_EMAIL:-$(git config user.email || echo methena@users.noreply.github.com)}" \
    commit -q -m "Set GitHub location to ${OWNER}/${NAME}" || true
if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  gh repo create "${OWNER}/${NAME}" $VIS --source . --remote origin --push \
     --description "Reproducible DNA-methylation pipeline: arrays & bisulfite-seq -> QC -> differential methylation, epigenetic clocks, regulation, Shiny"
  echo "Done: https://github.com/${OWNER}/${NAME}"
else
  echo
  echo "1) Create an EMPTY repository at https://github.com/new named '${NAME}' (no README/licence)."
  echo "2) Then run:"
  echo "   git remote add origin https://github.com/${OWNER}/${NAME}.git"
  echo "   git push -u origin main"
fi
