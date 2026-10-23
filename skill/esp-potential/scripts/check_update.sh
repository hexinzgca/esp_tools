#!/usr/bin/env bash
# Read-only preflight: compare the local esp_tools checkout with origin.
# Usage: check_update.sh [/path/to/esp_tools]
set -euo pipefail

REPO=${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}

if [[ ! -d "$REPO/.git" ]]; then
  echo "ERROR: not a git checkout: $REPO" >&2
  exit 2
fi

BRANCH=$(git -C "$REPO" symbolic-ref --short HEAD 2>/dev/null || echo master)
LOCAL=$(git -C "$REPO" rev-parse HEAD)
DIRTY=""
if [[ -n "$(git -C "$REPO" status --porcelain)" ]]; then
  DIRTY=" (working tree dirty)"
fi

if ! git -C "$REPO" fetch --quiet origin; then
  echo "WARNING: cannot fetch origin; local $LOCAL$DIRTY"
  exit 0
fi

if ! REMOTE=$(git -C "$REPO" rev-parse "origin/$BRANCH" 2>/dev/null); then
  echo "WARNING: origin/$BRANCH not found; local $LOCAL$DIRTY"
  exit 0
fi

if [[ "$LOCAL" == "$REMOTE" ]]; then
  echo "esp_tools up to date: $LOCAL$DIRTY"
  exit 0
fi

AHEAD=$(git -C "$REPO" rev-list --count "$REMOTE..$LOCAL")
BEHIND=$(git -C "$REPO" rev-list --count "$LOCAL..$REMOTE")
echo "esp_tools update: local $LOCAL$DIRTY"
echo "                  origin/$BRANCH $REMOTE (ahead $AHEAD, behind $BEHIND)"
if [[ "$BEHIND" -gt 0 ]]; then
  echo "run: git -C \"$REPO\" pull --ff-only"
fi
