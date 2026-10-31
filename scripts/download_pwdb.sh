#!/usr/bin/env bash
# Download the Pulse Wave Database (PWDB, synthetic) from Zenodo record 3275625, verify (MD5), unpack the WFDB waves.
# Only the files the converter needs (~190 MB); the large .mat files are skipped.
set -euo pipefail
cd "$(dirname "$0")/.."  # project root, so it works from any folder
mkdir -p data/raw/pwdb && cd data/raw/pwdb

URL=https://zenodo.org/api/records/3275625/files
# -C - resumes: just rerun the script if a download breaks
for f in PWs_wfdb.zip pwdb_haemod_params.csv pwdb_model_configs.csv pwdb_model_variations.csv; do
    curl -fL --retry 5 -# -o "$f" -C - "$URL/$f/content"
done

# MD5 sums from the Zenodo record
md5sum -c - <<'SUMS'
c4c3c3ba8163ee1f4f1495d7a4e70d73  PWs_wfdb.zip
43e1244665e6cee6b77501102404b70a  pwdb_haemod_params.csv
8c3b4f2f86386b72250766aedea9db52  pwdb_model_configs.csv
3f1987efbc131b2cab8e576807385d5b  pwdb_model_variations.csv
SUMS

unzip -q -o PWs_wfdb.zip -d PWs
# Once unpacked: rm PWs_wfdb.zip
