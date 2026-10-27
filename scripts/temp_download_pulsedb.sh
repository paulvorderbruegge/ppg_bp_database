#!/usr/bin/env bash
# Download PulseDB v2_0 Segment_Files (MIMIC + Vital) from Box, verify, unpack.
# Fedora: sudo dnf install 7zip   (older releases: p7zip p7zip-plugins). curl, awk, sha1sum are preinstalled.
set -euo pipefail
mkdir -p data/raw/pulsedb/Segment_Files && cd data/raw/pulsedb/Segment_Files

# -C - resumes: just rerun the script if a download breaks
dl() { curl -L -o "$1" -C - "$2"; }
dl PulseDB_MIMIC.zip.001 https://rutgers.box.com/shared/static/7l8n3tn9tr0602tdss1x7e3uliahlibp.001
dl PulseDB_MIMIC.zip.002 https://rutgers.box.com/shared/static/zco48rvz5dog72970679foen6hct15c8.002
dl PulseDB_MIMIC.zip.003 https://rutgers.box.com/shared/static/x22qpmelx6sz3wgkm5qyc0eis429361f.003
dl PulseDB_MIMIC.zip.004 https://rutgers.box.com/shared/static/xj25sqnluiz6s4z8tzzm5phk00ohp6e8.004
dl PulseDB_MIMIC.zip.005 https://rutgers.box.com/shared/static/dxus2lsoop02chaspnwipwrf0g4wmenr.005
dl PulseDB_MIMIC.zip.006 https://rutgers.box.com/shared/static/rts6sj441laenm2sy1qcemg7ke4om3j6.006
dl PulseDB_MIMIC.zip.007 https://rutgers.box.com/shared/static/vor4hjllld7a0c3nzef8uptbb4ut3koo.007
dl PulseDB_MIMIC.zip.008 https://rutgers.box.com/shared/static/a2qg2p4ebyrooji3z88djlokji65tlf3.008
dl PulseDB_MIMIC.zip.009 https://rutgers.box.com/shared/static/uh6kbiuqgnib5wakiv6o35gkpusyamc7.009
dl PulseDB_MIMIC.zip.010 https://rutgers.box.com/shared/static/h6eyhkkx48pf3ce3th1clwj43hn98j5c.010
dl PulseDB_MIMIC.zip.011 https://rutgers.box.com/shared/static/e93dp94hxpkas45yc59n289s2wvkafgi.011
dl PulseDB_MIMIC.zip.012 https://rutgers.box.com/shared/static/iuvyuw7dmlxvbjvt53dj49wqn3gelqni.012
dl PulseDB_MIMIC.zip.013 https://rutgers.box.com/shared/static/qxx6tjz8c3778601ib3icu6o1rranmc7.013
dl PulseDB_MIMIC.zip.014 https://rutgers.box.com/shared/static/ip2ninwqj8437l9fyffjprnk90ptnx9k.014
dl PulseDB_MIMIC.zip.015 https://rutgers.box.com/shared/static/yrtbo0lg8mjhaw624iw9bbhk1obbocwd.015
dl PulseDB_MIMIC.zip.016 https://rutgers.box.com/shared/static/wmzndowgfa5xi3tvtqahxkld3ngdyjds.016
dl PulseDB_Vital.zip.001 https://rutgers.box.com/shared/static/vtxoksmn7emeaxypb2prywgwscuefoqa.001
dl PulseDB_Vital.zip.002 https://rutgers.box.com/shared/static/euzkek7c3xoy62jisheuxqar7z5y8xig.002
dl PulseDB_Vital.zip.003 https://rutgers.box.com/shared/static/49lngo0benxfjw193jnqz9tctlyb3qam.003
dl PulseDB_Vital.zip.004 https://rutgers.box.com/shared/static/jf4fwgkmhry20mf5tcg9t0wxvky64um0.004
dl PulseDB_Vital.zip.005 https://rutgers.box.com/shared/static/2lgxysbskfuapsaan4jypvmm8316fdkc.005
dl PulseDB_Vital.zip.006 https://rutgers.box.com/shared/static/x27ktb4qsx43razwo4tjmxq9v1ro0x3y.006
dl PulseDB_Vital.zip.007 https://rutgers.box.com/shared/static/q0t36fikgf3pimhvnerwwnovfr0umtp8.007
dl PulseDB_Vital.zip.008 https://rutgers.box.com/shared/static/ihckx2g0f981g5yz2x8v5rgwndl6yebw.008
dl PulseDB_Vital.zip.009 https://rutgers.box.com/shared/static/y8j14h8tvi5b3du8nap9dnura1omfrk6.009
dl PulseDB_Vital.zip.010 https://rutgers.box.com/shared/static/fu0m9tx33jkxywq32shh0g8dg3not15u.010

# SHA-1 checks (from the README), "name: hash" -> "hash  name"
awk -F': ' '{print $2"  "$1}' <<'SUMS' | sha1sum -c -
PulseDB_MIMIC.zip.001: f3adc384962136eb93d5e12c73d1e2c742387df3
PulseDB_MIMIC.zip.002: b8db02f3e490c94e8b8b1e5c50bd127c0c340a65
PulseDB_MIMIC.zip.003: b32a27d7be4c8919b72a85b0bdafd716b0243f14
PulseDB_MIMIC.zip.004: 8e97e3811c42c6d5311d2a5151652be2b5b087d7
PulseDB_MIMIC.zip.005: 947ac26f686f3068cb726b3d158fb3469c203f68
PulseDB_MIMIC.zip.006: 8b6ab2773a3c7135d8cd71c9c33d7a480a1531e5
PulseDB_MIMIC.zip.007: 1648b9c50cf4b4949582955fa517c0a036982e03
PulseDB_MIMIC.zip.008: 4c9534c904d71cefafc6e8892d0eaeb6f6990113
PulseDB_MIMIC.zip.009: b7ad32b67abebec81861253ee1e1efac66a9527f
PulseDB_MIMIC.zip.010: 29fcc37cd04a9d099cb3a7566b07a3175b36e87f
PulseDB_MIMIC.zip.011: 70e585b97f2dc72130ea1c34ff0ba936a86e8e2c
PulseDB_MIMIC.zip.012: 5e4f38cc64cdb7938b144664f4863086bdf547f8
PulseDB_MIMIC.zip.013: 6acc91136ac53e232b4c83a6b3363b3077af2407
PulseDB_MIMIC.zip.014: 67b20543d8240ba6db7cf696638edf21c76ab26e
PulseDB_MIMIC.zip.015: 7e3edd06e61365b481233c60890309b43d57be29
PulseDB_MIMIC.zip.016: 837559cc18349e610d7950394e11be2d86d559c3
PulseDB_Vital.zip.001: 3e25f5f89e77b5619f911376f714facd1d14b95e
PulseDB_Vital.zip.002: 1ba93c9c4f1189f940be89db513b306ddbe93ffa
PulseDB_Vital.zip.003: 612ab2056288183bb2bba9724d6a66c18b15e71c
PulseDB_Vital.zip.004: f10ee6fe9f36b292ed5b313e6ea57357aa596c76
PulseDB_Vital.zip.005: 9a7c332ae6954a0e613b38b3ee522ab109e84736
PulseDB_Vital.zip.006: 8d7104a83e797d5c7026b0b1492d0f91fda6b426
PulseDB_Vital.zip.007: d60553c8cb46657f00f37c31494c5cde4b170a2b
PulseDB_Vital.zip.008: a64d6c051fe3f2bc0f6b438d6b9bcad7150e9c6f
PulseDB_Vital.zip.009: d98fb31b44364875c379fa82293ca14845b9f1c2
PulseDB_Vital.zip.010: 45c8a5bd810b1b5e70a63c3b9300490737fcabf4
SUMS

# The binary is called 7zz, 7z or 7za depending on the package
SEVENZ=$(command -v 7zz || command -v 7z || command -v 7za) || { echo "7-Zip not found"; exit 1; }

# 7-Zip picks up the remaining parts (.002, .003, ...) automatically
"$SEVENZ" x PulseDB_MIMIC.zip.001
"$SEVENZ" x PulseDB_Vital.zip.001

# Archives and unpacked files together need twice the space. Once everything is unpacked: rm PulseDB_*.zip.*
