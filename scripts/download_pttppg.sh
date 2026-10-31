#!/usr/bin/env bash
# Download PTT-PPG v1.1.0 (PhysioNet, open access, ~2.9 GB) to data/raw/pttppg and verify the SHA-256 sums.
set -euo pipefail
cd "$(dirname "$0")/.."  # project root, so it works from any folder
scripts/physionet_download.sh pulse-transit-time-ppg 1.1.0 data/raw/pttppg
