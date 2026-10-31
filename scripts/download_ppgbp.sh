#!/usr/bin/env bash
# Download the PPG-BP database (Liang et al. 2018, figshare 10.6084/m9.figshare.5459299.v5, ~2 MB), verify (MD5),
# unpack to data/raw/ppgbp/PPG-BP Database/Data File/.
set -euo pipefail
cd "$(dirname "$0")/.."  # project root, so it works from any folder
mkdir -p data/raw/ppgbp && cd data/raw/ppgbp

curl -fL --retry 5 -# -o "PPG-BP Database.zip" https://ndownloader.figshare.com/files/9441097
curl -fL --retry 5 -# -o "Table 1.xlsx" https://ndownloader.figshare.com/files/10395579

# MD5 sums from the figshare API
md5sum -c - <<'SUMS'
3b8eae44f45799aeb3c5f55af8589828  PPG-BP Database.zip
1d06355eb084a0c4e4506efe343bc3ff  Table 1.xlsx
SUMS

# The zip contains "Data File/..." at its top level
unzip -q -o "PPG-BP Database.zip" -d "PPG-BP Database"
rm "PPG-BP Database.zip"
