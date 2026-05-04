#!/bin/bash
#SBATCH --job-name=r_job
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=3-00:00:00
#SBATCH --mem=120G

set -euo pipefail

rscript=""
jobid=""
config_file=""

while getopts "r:j:c:" opt; do
  case $opt in
    r) rscript="$OPTARG" ;;
    j) jobid="$OPTARG" ;;
    c) config_file="$OPTARG" ;;
    *) echo "Usage: sbatch hpc/submit_r_job.sh -r <Rscript_file> -c <config_file> -j <jobid>" >&2; exit 1 ;;
  esac
done

if [[ -z "$rscript" || -z "$config_file" || -z "$jobid" ]]; then
  echo "Usage: sbatch hpc/submit_r_job.sh -r <Rscript_file> -c <config_file> -j <jobid>" >&2
  exit 1
fi
# example: sbatch hpc/submit_r_job.sh -r code/04_build_single_models.R -c config/config_build_single_models.csv -j GBM001

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"

mkdir -p logs

if [[ -n "$jobid" ]]; then
  scontrol update JobID="$SLURM_JOB_ID" JobName="$jobid"
fi

echo "Project directory: $PROJECT_DIR"
echo "R script: $rscript"
echo "Config file: $config_file"
echo "Job ID: $jobid"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "Start time: $(date)"

module load languages/R/4.3.3

if [[ ! -f "$rscript" ]]; then
  echo "Error: R script '$rscript' does not exist." >&2
  exit 1
fi

if [[ ! -f "$config_file" ]]; then
  echo "Error: config file '$config_file' does not exist." >&2
  exit 1
fi

Rscript "$rscript" "$config_file" "$jobid"

echo "End time: $(date)"