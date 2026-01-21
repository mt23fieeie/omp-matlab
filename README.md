# omp-matlab

Minimal MATLAB implementation of Orthogonal Matching Pursuit (OMP).

## Structure
- `src/omp.m`: OMP implementation (greedy atom selection + least-squares update)
- `demo/demo_sparse_recovery_sweep.m`: reproducible demo (single run + noise sweep)

## How to run
1. Open MATLAB at the repo root folder.
2. Make sure `src/` is on the MATLAB path (or `addpath("src")`).
3. Run:
   - `demo/demo_sparse_recovery_sweep.m`

## Output
The demo creates a `results/` folder and saves:
- residual plots: `results/residual_sigma*.png`
- logs: `results/run_sigma*.mat`

## Reproducibility
The demo uses `rng(0)` to keep the dictionary and sparse ground truth identical across noise levels.
