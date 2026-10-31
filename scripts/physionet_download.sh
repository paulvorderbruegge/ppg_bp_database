#!/usr/bin/env bash
# Download one open PhysioNet dataset file by file and check every file against PhysioNet's SHA256SUMS.txt.
# Used by download_pttppg.sh and download_bpgraphene.sh; can also be run directly:
#     scripts/physionet_download.sh <slug> <version> <target dir>
# Rerunning skips files that are already complete, so an interrupted download can simply be restarted.
set -euo pipefail
[ $# -eq 3 ] || { echo "usage: $0 <slug> <version> <target dir>"; exit 1; }
slug=$1 version=$2 target=$3
url="https://physionet.org/files/$slug/$version"

mkdir -p "$target" && cd "$target"
curl -fsSL -o SHA256SUMS.txt "$url/SHA256SUMS.txt"
echo "$slug $version: $(wc -l < SHA256SUMS.txt) files -> $PWD"

# Download every file whose checksum does not match yet, 4 at a time
export url
awk '{print $1, $2}' SHA256SUMS.txt | xargs -P 4 -n 2 sh -c '
    [ -f "$1" ] && echo "$0  $1" | sha256sum -c --status - && exit 0
    curl -fsSL --create-dirs --retry 5 -o "$1" "$url/$1" && echo "downloaded $1"
'

sha256sum -c --quiet SHA256SUMS.txt && echo "$slug: all checksums OK"
