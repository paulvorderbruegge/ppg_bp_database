"""Download a laptop-sized test set (~7.5 GB) that is representative of the full datasets, to data/raw/<dataset>.

    uv run python scripts/download_testset.py                 # all datasets
    uv run python scripts/download_testset.py pulsedb         # only these
    uv run python scripts/download_testset.py --pulsedb-gb 3  # smaller PulseDB sample

What is downloaded:
- ppgbp, pwdb, mimicbp: complete (download_<dataset>.sh), archives deleted after unpacking.
- pttppg, bpgraphene: complete except files the converters never read (pttppg: the CSV copies of the WFDB records in
  csv/, except subjects_info.csv; bpgraphene: the *_bioz.csv bioimpedance). Every file is checked against PhysioNet's
  SHA-256 sum; SHA256SUMS.txt is trimmed to the downloaded files, so `sha256sum -c SHA256SUMS.txt` still works.
- pulsedb: a random sample of whole subject files (seed --seed). The budget is split between the MIMIC and VitalDB
  parts in proportion to their unpacked size in the full dataset (MIMIC 459 GB : VitalDB 177 GB), so the mix of ICU
  and surgery data matches. Within a part, the sample is stratified by file size (= recording length, quintiles),
  with the same number of subjects per quintile, so short and long recordings appear as often as in the full
  dataset (see _stratified). Each file is streamed out of the split zip with HTTP range requests and
  checked against the CRC-32 stored in the zip. The chosen files are listed in data/raw/pulsedb/TESTSET.

Each dataset folder gets a TESTSET marker file. A folder that holds data without a TESTSET or SAMPLE marker (i.e. a
full download) is never touched; a finished dataset is skipped on a rerun.
"""
import argparse
import hashlib
import json
import random
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile
import zlib
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import download_sample as sample  # noqa: E402  (RemoteFile, PULSEDB_PARTS)

import numpy as np  # noqa: E402

from ppg_bp_database.utils import RAW_DIR  # noqa: E402

SCRIPTS = Path(__file__).parent
MARKER = "TESTSET"
HEADERS = {"User-Agent": "curl/8"}
CHUNK = 4 * 2**20

# PhysioNet: (slug, version, regex of files NOT to download)
PHYSIONET = {
    "pttppg": ("pulse-transit-time-ppg", "1.1.0", r"^csv/(?!subjects_info\.csv$)"),
    "bpgraphene": ("bp-graphene-bioimpedance", "1.0.0", r"_bioz\.csv$"),
}
# unpacked size of the full PulseDB parts (GB), from the zip directory; used to split the budget
PULSEDB_FULL_GB = {"PulseDB_MIMIC": 459.0, "PulseDB_Vital": 177.2}
PULSEDB_STRATA = 5  # size quintiles, see _stratified


def _open(url: str, start: int | None = None, end: int | None = None):
    headers = dict(HEADERS)
    if start is not None:
        headers["Range"] = f"bytes={start}-{end}"
    return urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120)


def _physionet_file(url: str, name: str, sha256: str, out: Path) -> int:
    dest = out / name
    if dest.is_file() and hashlib.sha256(dest.read_bytes()).hexdigest() == sha256:
        return dest.stat().st_size
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(5):
        h = hashlib.sha256()
        try:
            with _open(f"{url}/{name}") as r, dest.open("wb") as f:
                while chunk := r.read(CHUNK):
                    h.update(chunk)
                    f.write(chunk)
        except OSError as e:
            if attempt == 4:
                raise RuntimeError(f"{name}: {e}") from e
            continue
        if h.hexdigest() == sha256:
            return dest.stat().st_size
    raise RuntimeError(f"{name}: SHA-256 mismatch after 5 attempts")


def physionet(name: str, out: Path) -> str:
    slug, version, exclude = PHYSIONET[name]
    url = f"https://physionet.org/files/{slug}/{version}"
    with _open(f"{url}/SHA256SUMS.txt") as r:
        lines = r.read().decode().splitlines()
    keep = [line for line in lines if line.strip() and not re.search(exclude, line.split(maxsplit=1)[1])]
    print(f"  {len(keep)} of {len(lines)} files ({len(lines) - len(keep)} excluded: {exclude})")
    with ThreadPoolExecutor(8) as pool:
        sizes = list(pool.map(lambda line: _physionet_file(url, line.split(maxsplit=1)[1], line.split()[0], out), keep))
    (out / "SHA256SUMS.txt").write_text("\n".join(keep) + "\n")
    return (f"PhysioNet {slug} {version}: {len(keep)} of {len(lines)} files ({sum(sizes) / 1e9:.2f} GB), "
            f"excluded (not read by the converter): {exclude}")


def full(name: str, out: Path) -> str:
    """Run download_<name>.sh (it writes to data/raw/<name>), then delete the archives it keeps."""
    subprocess.run([str(SCRIPTS / f"download_{name}.sh")], check=True)
    for z in out.glob("*.zip"):
        z.unlink()
    return f"complete (download_{name}.sh, archives deleted after unpacking)"


class _Parts:
    """Byte ranges of a split zip (several URLs = one concatenated file), streamed with one request per part."""

    def __init__(self, remote: sample.RemoteFile):
        self.parts = remote.parts  # (url, offset, size)

    def stream(self, start: int, length: int):
        end = start + length
        for url, offset, size in self.parts:
            a, b = max(start, offset), min(end, offset + size)
            if a < b:
                with _open(url, a - offset, b - offset - 1) as r:
                    if r.status != 206:
                        raise RuntimeError("server ignored the range request")
                    while chunk := r.read(CHUNK):
                        yield chunk


def _extract_streamed(remote: sample.RemoteFile, parts: _Parts, info, dest: Path) -> None:
    """Extract one zip member by streaming its compressed bytes in one go (deflate or stored), checking the CRC."""
    remote.seek(info.header_offset)
    header = remote.read(30)
    assert header[:4] == b"PK\x03\x04", f"{info.filename}: no local file header"
    data_start = info.header_offset + 30 + int.from_bytes(header[26:28], "little") + int.from_bytes(header[28:30], "little")
    if info.compress_type not in (0, 8):
        raise RuntimeError(f"{info.filename}: unsupported compression {info.compress_type}")
    inflate = zlib.decompressobj(-15) if info.compress_type == 8 else None
    crc, n = 0, 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    with tmp.open("wb") as f:
        for chunk in parts.stream(data_start, info.compress_size):
            data = inflate.decompress(chunk) if inflate else chunk
            crc, n = zlib.crc32(data, crc), n + len(data)
            f.write(data)
        if inflate:
            data = inflate.flush()
            crc, n = zlib.crc32(data, crc), n + len(data)
            f.write(data)
    if crc != info.CRC or n != info.file_size:
        tmp.unlink()
        raise RuntimeError(f"{info.filename}: CRC / size mismatch")
    tmp.rename(dest)


def _stratified(infos: list, budget: float, seed: int) -> list:
    """Subject files sorted by size, cut into PULSEDB_STRATA strata of equal count (recording-length quantiles), and
    the same number k of random files taken from each stratum, k as large as the budget allows. The sample then has
    the same distribution of recording lengths as the full dataset; a plain random draw within a budget does not
    (a few very long recordings use up the budget, or filling the rest with small files skews it to short ones)."""
    by_size = sorted(infos, key=lambda i: (i.file_size, i.filename))
    strata = []
    for s, chunk in enumerate(np.array_split(np.arange(len(by_size)), PULSEDB_STRATA)):
        chunk = [by_size[j] for j in chunk]
        random.Random(f"{seed}-{s}").shuffle(chunk)
        strata.append(chunk)
    k = 0
    while k < len(strata[0]) and sum(i.file_size for st in strata for i in st[:k + 1]) <= budget:
        k += 1
    return sorted((i for st in strata for i in st[:k]), key=lambda i: i.filename)


def pulsedb(out: Path, budget_gb: float, seed: int) -> str:
    sample.MAX_TOTAL_BYTES = 2**40  # the limit is for the sample script; here only the zip directories are read
    total_full = sum(PULSEDB_FULL_GB.values())
    chosen, lines = {}, []
    for part, ids in sample.PULSEDB_PARTS.items():
        remote = sample.RemoteFile([f"https://rutgers.box.com/shared/static/{p}" for p in ids])
        infos = sorted((i for i in zipfile.ZipFile(remote).infolist() if i.filename.endswith(".mat")),
                       key=lambda i: i.filename)
        budget = budget_gb * 1e9 * PULSEDB_FULL_GB[part] / total_full
        picked = _stratified(infos, budget, seed)
        used = sum(i.file_size for i in picked)
        chosen[part] = [Path(i.filename).name for i in picked]
        print(f"  {part}: {len(picked)} of {len(infos)} subjects, {used / 1e9:.2f} GB unpacked (budget {budget / 1e9:.2f})")
        # files of an earlier selection (other seed / budget) are removed
        names = {Path(i.filename).name for i in picked}
        for old in (out / "Segment_Files" / part).glob("*.mat"):
            if old.name not in names:
                old.unlink()
        parts = _Parts(remote)
        for k, info in enumerate(picked, 1):
            dest = out / "Segment_Files" / info.filename
            if dest.is_file() and dest.stat().st_size == info.file_size:
                continue
            print(f"    [{k}/{len(picked)}] {info.filename} ({info.file_size / 1e6:.0f} MB)", flush=True)
            _extract_streamed(remote, parts, info, dest)
        lines.append(f"{part}: {len(picked)} of {len(infos)} subject files, {used / 1e9:.2f} of "
                     f"{sum(i.file_size for i in infos) / 1e9:.1f} GB unpacked")
    (out / "subjects.json").write_text(json.dumps({"seed": seed, "budget_gb": budget_gb, "files": chosen}, indent=1))
    return (f"random sample of subject files stratified by size quintile, seed {seed}, budget {budget_gb:g} GB "
            f"split by unpacked size "
            f"(list: subjects.json)\n" + "\n".join(lines))


DATASETS = ("ppgbp", "pttppg", "bpgraphene", "mimicbp", "pwdb", "pulsedb")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="*", choices=[[], *DATASETS], metavar="dataset")
    parser.add_argument("--pulsedb-gb", type=float, default=5.0, help="unpacked size of the PulseDB sample (GB)")
    parser.add_argument("--seed", type=int, default=0, help="seed of the PulseDB subject sample")
    args = parser.parse_args()

    failed = []
    for name in args.datasets or DATASETS:
        print(f"--- {name}", flush=True)
        out = RAW_DIR / name
        if (out / MARKER).exists():
            print("  already complete (delete the folder to download again)")
            continue
        if out.is_dir() and any(out.iterdir()):
            if not (out / sample.MARKER).exists() and not (out / f"{MARKER}.partial").exists():
                print(f"  SKIPPED: {out} holds data that is neither a sample nor a test set (a full download?)")
                failed.append(name)
                continue
            if (out / sample.MARKER).exists():
                shutil.rmtree(out)
        out.mkdir(parents=True, exist_ok=True)
        (out / f"{MARKER}.partial").write_text("download_testset.py did not finish this dataset yet; rerun it\n")
        try:
            if name in PHYSIONET:
                how = physionet(name, out)
            elif name == "pulsedb":
                how = pulsedb(out, args.pulsedb_gb, args.seed)
            else:
                how = full(name, out)
        except Exception as e:  # keep going with the other datasets; a rerun resumes
            print(f"  FAILED: {type(e).__name__}: {e}")
            failed.append(name)
            continue
        (out / MARKER).write_text(f"Test set written by scripts/download_testset.py on {date.today()}.\n{how}\n")
        (out / f"{MARKER}.partial").unlink()
        print(f"  {how}")
    if failed:
        sys.exit(f"failed: {', '.join(failed)} (rerun to resume)")


if __name__ == "__main__":
    main()
