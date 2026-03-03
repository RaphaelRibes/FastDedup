#!/bin/bash
#SBATCH --job-name=bench_fdedup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal

pixi run cargo build --release
cp ../target/release/fdedup /bin/
pixi run -e bench ./benchmark.sh
pixi run -e bench python plotit.py