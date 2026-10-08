#!/bin/bash
#SBATCH --job-name=figshare-plots
#SBATCH --account=torch_pr_126_general
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --chdir=/scratch/ycc520/flyem
#SBATCH --output=/scratch/ycc520/flyem/int/figshare_pooled_upload/slurm-%j.log
#SBATCH --open-mode=append

set -euo pipefail
for attempt in 1 2 3 4; do
  result=0
  /usr/bin/python -u src/upload_figshare_plots.py \
    --token-file /scratch/ycc520/flyem/int/FS_TOKEN \
    --staging /scratch/ycc520/flyem/int/figshare_pooled_upload --workers 4 || result=$?
  if [ "$result" -eq 0 ]; then exit 0; fi
  if [ "$result" -ne 75 ] || [ "$attempt" -eq 4 ]; then exit "$result"; fi
  echo "Figshare verification delayed; resuming from checkpoints in 60 seconds."
  sleep 60
done
