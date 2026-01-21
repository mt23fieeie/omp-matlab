# `omp-matlab`

A minimal MATLAB implementation of **Orthogonal Matching Pursuit (OMP)** for sparse recovery.

Given a dictionary $A \in \mathbb{R}^{n \times m}$ and an observation $b \in \mathbb{R}^n$, recover a **$K$-sparse** coefficient vector $x$ such that

$$
b \approx A x.
$$

## Repository Structure

This repository contains:

- `src/omp.m` — Core OMP implementation (greedy atom selection + least-squares update)  
- `demo/demo_sparse_recovery_sweep.m` — Reproducible demo script (single run + noise sweep)

## Quick Start

Open MATLAB at the repository root and run:

```matlab
addpath("src")
run("demo/demo_sparse_recovery_sweep.m")
```

The demo will create a `results/` folder and save:

- Residual plots: `results/residual_sigma*.png`
- Logs: `results/run_sigma*.mat`

## Reproducibility

The demo uses:

```matlab
rng(0)
```

to ensure that:

- The dictionary $A$
- The ground-truth sparse vector $x_0$

are identical across different noise levels (fair comparison). Only the additive Gaussian noise changes.

## Example Results (Noise Sweep)

### Synthetic Setup

- Random Gaussian dictionary $(n = 64,\; m = 128)$
- Sparsity $K = 6$
- Observation model:

$$
b = A x_0 + \sigma \varepsilon, \qquad \varepsilon \sim \mathcal{N}(0, I).
$$

## Metrics

To evaluate the recovery performance, we use the following three indicators:

* **Hit Count**: The number of correctly recovered support indices, representing the intersection between the ground-truth support ($S_0$) and the estimated support ($\hat{S}$):
    $$|S_0 \cap \hat{S}|$$

* **Relative Coefficient Error (rel_err_x)**: Measures the distance between the estimated sparse vector $\hat{x}$ and the ground-truth $x_0$:
    $$\frac{\|\hat{x} - x_0\|_2}{\|x_0\|_2}$$

* **Relative Residual (rel_resid)**: Measures how well the recovered signal approximates the observation $b$ in the measurement space:
    $$\frac{\|b - A \hat{x}\|_2}{\|b\|_2}$$

## Results Table

| sigma | hit | rel_err_x | rel_resid |
|------:|---:|----------:|----------:|
| 0.01  | 6  | 1.84e-03  | 4.98e-03  |
| 0.1   | 6  | 1.84e-02  | 4.99e-02  |
| 1     | 5  | 2.66e-01  | 4.40e-01  |
| 10    | 0  | 4.75e+00  | 7.01e-01  |

## Observation

OMP:

- Recovers the full support at **low noise**
- Gradually degrades at **moderate noise**
- Fails at **high noise**, where correlation-based greedy selection becomes noise-dominated

## Notes

- `normalize_cols = true` (default): Normalizes dictionary columns to avoid selection bias caused by unequal norms.
- The demo script includes:
  - A **single-run mode**
  - A **noise-sweep mode** for controlled comparison

