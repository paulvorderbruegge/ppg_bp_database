"""Download a tiny sample of every dataset (a few records each) to data/raw/<dataset>, to test the pipeline end to
end. Not for analysis.

    uv run python scripts/download_sample.py              # all datasets
    uv run python scripts/download_sample.py pwdb mimicbp # only these
    uv run python scripts/structure.py
    uv run python scripts/harmonize.py

A sample folder is marked with a file named SAMPLE. The script only replaces a dataset folder that is empty or
has this marker, so a full download (download_<dataset>.sh) is never deleted. To replace a sample with the full
data: delete the folder, then run download_<dataset>.sh.

Only single files are fetched. Plain files take one request. For files inside a zip (also PulseDB's split zip),
HTTP range requests read only the zip's directory and the wanted members, never the whole archive.
Hard limits: a file (compressed size) above MAX_FILE_BYTES is refused, and the run stops before the total
download exceeds MAX_TOTAL_BYTES. No checksum checks, because only parts are downloaded and
pttppg's subjects_info.csv is trimmed to the downloaded records.
"""
import argparse
import csv
import io
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

from ppg_bp_database.utils import RAW_DIR as OUT_DIR

MARKER = "SAMPLE"  # marks a dataset folder as written by this script

MAX_FILE_BYTES = 15 * 2**20
MAX_TOTAL_BYTES = 100 * 2**20
BLOCK = 64 * 2**10  # granularity of range requests into a zip
HEADERS = {"User-Agent": "curl/8"}

downloaded = 0  # bytes, over the whole run


class TooLarge(Exception):
    pass


def _get(url: str, start: int | None = None, end: int | None = None) -> tuple[bytes, str, int | None]:
    """GET (bytes start..end inclusive, or the whole file). Returns (body, final URL, total size from Content-Range)."""
    global downloaded
    headers = dict(HEADERS)
    if start is not None:
        headers["Range"] = f"bytes={start}-{end}"
        if downloaded + (end - start + 1) > MAX_TOTAL_BYTES:
            raise TooLarge(f"total download would exceed {MAX_TOTAL_BYTES / 2**20:.0f} MB")
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60) as r:
        if start is not None and r.status != 206:
            raise RuntimeError(f"{url}: server ignored the range request (status {r.status})")
        total = int(r.headers["Content-Range"].rsplit("/", 1)[1]) if start is not None else None
        body = r.read(MAX_TOTAL_BYTES - downloaded + 1)
        downloaded += len(body)
        if downloaded > MAX_TOTAL_BYTES:
            raise TooLarge(f"total download exceeded {MAX_TOTAL_BYTES / 2**20:.0f} MB")
        return body, r.url, total


def size_of(url: str) -> tuple[int, str]:
    """(size in bytes, final URL after redirects), with a 1-byte range request."""
    _, final, total = _get(url, 0, 0)
    return total, final


def fetch(url: str, dest: Path) -> None:
    """Download one whole file, refused if larger than MAX_FILE_BYTES."""
    size, final = size_of(url)
    if size > MAX_FILE_BYTES:
        raise TooLarge(f"{dest.name}: {size / 2**20:.1f} MB > {MAX_FILE_BYTES / 2**20:.0f} MB")
    body, _, _ = _get(final, 0, size - 1)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(body)
    print(f"  {dest.relative_to(OUT_DIR)}  ({size / 2**10:.0f} kB)")


class RemoteFile(io.RawIOBase):
    """Read-only, seekable view of a remote file (or of the concatenation of several parts, e.g. a split zip),
    read with HTTP range requests in blocks of BLOCK bytes. Enough for zipfile.ZipFile."""

    def __init__(self, urls: list[str]):
        self.parts = []  # (resolved URL, offset, size)
        offset = 0
        for url in urls:
            size, final = size_of(url)
            self.parts.append((final, offset, size))
            offset += size
        self.size, self.pos, self.cache = offset, 0, {}

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, offset, whence=io.SEEK_SET):
        self.pos = {io.SEEK_SET: 0, io.SEEK_CUR: self.pos, io.SEEK_END: self.size}[whence] + offset
        return self.pos

    def _block(self, i: int) -> bytes:
        if i not in self.cache:
            start, end = i * BLOCK, min((i + 1) * BLOCK, self.size)
            chunks = []
            for url, offset, size in self.parts:
                a, b = max(start, offset), min(end, offset + size)
                if a < b:
                    chunks.append(_get(url, a - offset, b - offset - 1)[0])
            self.cache[i] = b"".join(chunks)
        return self.cache[i]

    def readinto(self, buf):
        n = min(len(buf), max(0, self.size - self.pos))
        out = bytearray()
        while len(out) < n:
            block = self._block((self.pos + len(out)) // BLOCK)
            start = (self.pos + len(out)) % BLOCK
            out += block[start:start + n - len(out)]
        buf[:n] = out
        self.pos += n
        return n


def extract(z: zipfile.ZipFile, name: str, dest_dir: Path) -> None:
    """Extract one zip member to dest_dir/<its path in the zip>, refused if larger than MAX_FILE_BYTES."""
    info = z.getinfo(name)
    if info.compress_size > MAX_FILE_BYTES:
        raise TooLarge(f"{name}: {info.compress_size / 2**20:.1f} MB > {MAX_FILE_BYTES / 2**20:.0f} MB")
    path = z.extract(info, dest_dir)
    print(f"  {Path(path).relative_to(OUT_DIR)}  ({info.compress_size / 2**10:.0f} kB compressed)")


def remote_zip(*urls: str) -> zipfile.ZipFile:
    return zipfile.ZipFile(RemoteFile(list(urls)))


# --- one function per dataset; each writes the minimum its converter needs ------------------------------------------

def ppgbp(out: Path) -> None:
    """Excel table + the 3 segments of 2 subjects."""
    z = remote_zip("https://ndownloader.figshare.com/files/9441097")
    names = z.namelist()
    extract(z, next(n for n in names if n.endswith("PPG-BP dataset.xlsx")), out / "PPG-BP Database")
    txt = sorted(n for n in names if "/0_subject/" in n and n.endswith(".txt"))
    subjects = sorted({Path(n).stem.split("_")[0] for n in txt}, key=int)[:2]
    for n in txt:
        if Path(n).stem.split("_")[0] in subjects:
            extract(z, n, out / "PPG-BP Database")


def pttppg(out: Path) -> None:
    """Two records (s1 sitting and walking); subjects_info.csv trimmed to them, because the converter reads every
    record listed there."""
    url = "https://physionet.org/files/pulse-transit-time-ppg/1.1.0"
    records = ["s1_sit", "s1_walk"]
    for rec in records:
        for ext in ("hea", "dat"):
            fetch(f"{url}/{rec}.{ext}", out / f"{rec}.{ext}")
    info = out / "csv/subjects_info.csv"
    fetch(f"{url}/csv/subjects_info.csv", info)
    with info.open(newline="") as f:
        rows = list(csv.reader(f))
    col = rows[0].index("record")
    with info.open("w", newline="") as f:
        csv.writer(f).writerows([rows[0]] + [r for r in rows[1:] if r[col] in records])


def bpgraphene(out: Path) -> None:
    """One trial (subject 1, day 1, baseline, trial 1), without the unused bioimpedance."""
    url = "https://physionet.org/files/bp-graphene-bioimpedance/1.0.0"
    trial = "subject1_day1/setup01_baseline/data_trial01"
    for signal in ("ppg", "finapresPPG", "finapresBP"):
        fetch(f"{url}/{trial}_{signal}.csv", out / f"{trial}_{signal}.csv")


def mimicbp(out: Path) -> None:
    """ppg, abp and labels of 2 subjects, each from its own zip."""
    url = "https://dataverse.harvard.edu/api/access/datafile/{}"
    zips = {"ppg": 7574732, "abp": 7574731, "labels": 7574727}
    ppg = remote_zip(url.format(zips["ppg"]))
    names = sorted(n for n in ppg.namelist() if n.endswith("_ppg.npy"))[:2]
    subjects = [Path(n).name.removesuffix("_ppg.npy") for n in names]
    for kind, file_id in zips.items():
        z = ppg if kind == "ppg" else remote_zip(url.format(file_id))
        for subject in subjects:
            member = next(n for n in z.namelist() if Path(n).name == f"{subject}_{kind}.npy")
            extract(z, member, out / kind)  # same folders as download_mimicbp.sh


def pwdb(out: Path) -> None:
    """The three small subject tables (all subjects) + the WFDB records of 3 virtual subjects."""
    url = "https://zenodo.org/api/records/3275625/files/{}/content"
    for name in ("pwdb_haemod_params.csv", "pwdb_model_configs.csv", "pwdb_model_variations.csv"):
        fetch(url.format(name), out / name)
    z = remote_zip(url.format("PWs_wfdb.zip"))
    headers = sorted(n for n in z.namelist() if Path(n).name.startswith("pwdb") and n.endswith(".hea"))[:3]
    for h in headers:
        extract(z, h, out / "PWs")
        extract(z, h.removesuffix(".hea") + ".dat", out / "PWs")


PULSEDB_PARTS = {  # same links as download_pulsedb.sh
    "PulseDB_MIMIC": ["7l8n3tn9tr0602tdss1x7e3uliahlibp.001", "zco48rvz5dog72970679foen6hct15c8.002",
                      "x22qpmelx6sz3wgkm5qyc0eis429361f.003", "xj25sqnluiz6s4z8tzzm5phk00ohp6e8.004",
                      "dxus2lsoop02chaspnwipwrf0g4wmenr.005", "rts6sj441laenm2sy1qcemg7ke4om3j6.006",
                      "vor4hjllld7a0c3nzef8uptbb4ut3koo.007", "a2qg2p4ebyrooji3z88djlokji65tlf3.008",
                      "uh6kbiuqgnib5wakiv6o35gkpusyamc7.009", "h6eyhkkx48pf3ce3th1clwj43hn98j5c.010",
                      "e93dp94hxpkas45yc59n289s2wvkafgi.011", "iuvyuw7dmlxvbjvt53dj49wqn3gelqni.012",
                      "qxx6tjz8c3778601ib3icu6o1rranmc7.013", "ip2ninwqj8437l9fyffjprnk90ptnx9k.014",
                      "yrtbo0lg8mjhaw624iw9bbhk1obbocwd.015", "wmzndowgfa5xi3tvtqahxkld3ngdyjds.016"],
    "PulseDB_Vital": ["vtxoksmn7emeaxypb2prywgwscuefoqa.001", "euzkek7c3xoy62jisheuxqar7z5y8xig.002",
                      "49lngo0benxfjw193jnqz9tctlyb3qam.003", "jf4fwgkmhry20mf5tcg9t0wxvky64um0.004",
                      "2lgxysbskfuapsaan4jypvmm8316fdkc.005", "x27ktb4qsx43razwo4tjmxq9v1ro0x3y.006",
                      "q0t36fikgf3pimhvnerwwnovfr0umtp8.007", "ihckx2g0f981g5yz2x8v5rgwndl6yebw.008",
                      "y8j14h8tvi5b3du8nap9dnura1omfrk6.009", "fu0m9tx33jkxywq32shh0g8dg3not15u.010"],
}


def pulsedb(out: Path) -> None:
    """The smallest subject file of each part (MIMIC, VitalDB). The split archive (.zip.001, ...) is one zip cut
    into byte chunks, so the parts are read as one concatenated remote file."""
    for parts in PULSEDB_PARTS.values():
        z = remote_zip(*(f"https://rutgers.box.com/shared/static/{p}" for p in parts))
        mats = [i for i in z.infolist() if i.filename.endswith(".mat")]
        extract(z, min(mats, key=lambda i: i.compress_size).filename, out / "Segment_Files")


DATASETS = {f.__name__: f for f in (ppgbp, pttppg, bpgraphene, mimicbp, pwdb, pulsedb)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("datasets", nargs="*", choices=[[], *DATASETS], metavar="dataset")
    names = parser.parse_args().datasets or list(DATASETS)

    failed = []
    for name in names:
        print(f"--- {name}")
        out = OUT_DIR / name
        if out.is_dir() and any(out.iterdir()) and not (out / MARKER).exists():
            print(f"SKIPPED: {out} already holds data that is not a sample (delete it first to get a sample)")
            failed.append(name)
            continue
        shutil.rmtree(out, ignore_errors=True)  # a sample is always downloaded fresh
        out.mkdir(parents=True)
        try:
            DATASETS[name](out)
            (out / MARKER).write_text("Small test sample written by scripts/download_sample.py, not the full dataset.\n")
        except Exception as e:
            print(f"SKIPPED: {type(e).__name__}: {e}")
            shutil.rmtree(out)
            failed.append(name)
            if isinstance(e, TooLarge) and downloaded >= MAX_TOTAL_BYTES:
                break
    print(f"downloaded {downloaded / 2**20:.1f} MB in total")
    if failed:
        sys.exit(f"skipped: {', '.join(failed)}")


if __name__ == "__main__":
    main()
