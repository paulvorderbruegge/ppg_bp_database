#!/usr/bin/env bash
# Download the parts of the MIMIC-III Clinical Database v1.4 needed for demographics (age, sex, height, weight) of
# MIMIC-BP and the MIMIC part of PulseDB to data/raw/mimic3, verified with PhysioNet's SHA-256 checksums.
#     PHYSIONET_USER=<username> scripts/download_mimic3.sh      (or: scripts/download_mimic3.sh <username>)
# MIMIC-III is credentialed: your PhysioNet account must have signed the data use agreement. The password is asked
# once (or taken from ~/.netrc if that has an entry for physionet.org).
# Uses aria2c if installed (Fedora: sudo dnf install aria2): PhysioNet limits the speed per connection, so aria2c
# downloads each file over several connections at once (MIMIC3_CONNECTIONS, default 8, max 16). Else wget.
#
# Files (~4.0 GB, CHARTEVENTS alone 4.0 GB; it stays gzipped):
#     PATIENTS, ADMISSIONS          sex, date of birth, admission times (age at the time of the recording)
#     CHARTEVENTS, D_ITEMS          height / weight measurements (+ item labels, to print what was used)
#     mimic3wdb-matched/RECORDS-waveforms   open access: start time of every waveform record of the MIMIC-III
#                                   Waveform Database Matched Subset, to find the admission a recording comes from
# Rerunning skips complete files and resumes broken ones.
set -euo pipefail
cd "$(dirname "$0")/.."  # project root, so it works from any folder
mkdir -p data/raw/mimic3/mimic3wdb-matched && cd data/raw/mimic3

url=https://physionet.org/files/mimiciii/1.4
user=${1:-${PHYSIONET_USER:-}}
files=(SHA256SUMS.txt PATIENTS.csv.gz ADMISSIONS.csv.gz D_ITEMS.csv.gz CHARTEVENTS.csv.gz)

# aria2c / wget, not curl: PhysioNet answers curl's requests to credentialed files with 403 before checking the
# login. aria2c sends wget's user agent to be safe.
have_netrc=false; grep -qs "physionet.org" ~/.netrc && have_netrc=true
if ! $have_netrc && [ -z "$user" ]; then
    echo "usage: PHYSIONET_USER=<username> $0   (or: $0 <username>)"; exit 1
fi
fail_msg="download failed (401/403: wrong username/password, or MIMIC-III data use agreement not signed)"

if command -v aria2c > /dev/null; then
    n=${MIMIC3_CONNECTIONS:-8}
    # The URLs (and the password) go in through stdin, not the command line, so the password is not visible in ps.
    input=""
    if ! $have_netrc; then
        read -rsp "Password for $user@physionet.org: " pw; echo
    fi
    for f in "${files[@]}"; do
        input+="$url/$f"$'\n'
        $have_netrc || input+="  http-user=$user"$'\n'"  http-passwd=$pw"$'\n'
    done
    # -c resumes broken downloads (also ones started by wget) and skips complete files; -j: files in parallel,
    # -x/-s: connections per file. --auto-file-renaming=false: never write CHARTEVENTS.1.csv.gz next to an old file.
    aria2c -i - -c -j 4 -x "$n" -s "$n" -k 1M --max-tries=5 --retry-wait=5 --file-allocation=none \
        --auto-file-renaming=false --allow-overwrite=false --user-agent="Wget/1.21.4" \
        --summary-interval=30 --console-log-level=warn <<< "$input" \
        || { echo "$fail_msg"; exit 1; }
    unset pw input
else
    command -v wget > /dev/null || { echo "neither aria2c nor wget found (Fedora: sudo dnf install aria2)"; exit 1; }
    echo "aria2c not found, using wget (one connection per file, slow; Fedora: sudo dnf install aria2)"
    auth=(); $have_netrc || auth=(--user "$user" --ask-password)  # wget reads ~/.netrc itself
    # One wget call for all files, so the password is asked only once. -c resumes broken downloads and skips
    # complete files.
    wget -c --tries=5 "${auth[@]}" "${files[@]/#/$url/}" || { echo "$fail_msg"; exit 1; }
fi

for f in "${files[@]:1}"; do
    sum=$(awk -v f="$f" '$2 == f {print $1}' SHA256SUMS.txt)
    [ -n "$sum" ] || { echo "$f is not listed in SHA256SUMS.txt"; exit 1; }
    echo "$sum  $f" | sha256sum -c - || { echo "delete $f and rerun"; exit 1; }
done

# Open access (no login). Checksum taken from mimic3wdb-matched/1.0/SHA256SUMS.txt (172 MB, so not downloaded here).
f=mimic3wdb-matched/RECORDS-waveforms
sum=48de2a92e38d1ea76281352f26f0cf8ea5022c801ac136c7532648607cfbef9d
if ! { [ -f "$f" ] && echo "$sum  $f" | sha256sum -c --status -; }; then
    curl -fsSL --retry 5 -o "$f" "https://physionet.org/files/mimic3wdb-matched/1.0/RECORDS-waveforms"
fi
echo "$sum  $f" | sha256sum -c -

echo "MIMIC-III OK. The first structure run of mimicbp / pulsedb extracts height and weight from CHARTEVENTS once."
