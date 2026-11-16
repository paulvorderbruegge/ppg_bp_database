#!/bin/bash
# Full download + structure + harmonize of all datasets except PulseDB, on the LRZ Linux Cluster (CoolMUC-4).
# Submit from the project folder:  sbatch jobs/full_nopulsedb.sh
# serial_std is used because its compute nodes have internet access (needed for the downloads).
#SBATCH -J ppg_nopulsedb
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00

module load slurm_setup
set -euo pipefail
UV=$HOME/.local/bin/uv

scripts/download_all.sh ppgbp pttppg bpgraphene mimicbp pwdb
$UV run python scripts/structure.py ppgbp pttppg bpgraphene mimicbp pwdb
$UV run python scripts/harmonize.py
