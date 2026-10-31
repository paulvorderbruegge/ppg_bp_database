#!/usr/bin/env bash
# Download the graphene bioimpedance BP dataset v1.0.0 (PhysioNet, open access, ~2.3 GB) to data/raw/bpgraphene
# and verify the SHA-256 sums. All files are downloaded, including the bioimpedance (*_bioz.csv) that the converter
# does not use yet, so the raw folder is a complete copy of the dataset.
set -euo pipefail
cd "$(dirname "$0")/.."  # project root, so it works from any folder
scripts/physionet_download.sh bp-graphene-bioimpedance 1.0.0 data/raw/bpgraphene
