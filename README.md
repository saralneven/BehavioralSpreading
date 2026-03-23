# Neven2026_BehavioralSpreading

Code and data for:

**Neven, S.** (2026). *Behavioral spreading in reef fish groups: a social contagion model of startle responses.* University of Amsterdam.

This repository contains all code and processed data needed to reproduce the analyses, model fitting, simulations, and figures reported in the paper. The study models the spreading of startle responses through reef fish groups using a hazard-based social contagion framework in which a fish's response probability depends on the timing, proximity, and speed of its already-responding neighbours.

---

## Repository structure

```
data/
  all_observations.csv                    — all fish per video (output of 01_video_processing.ipynb)
  filtered_observations.csv              — social responders only (SRL, SRNL, NRL_r, NRNL_r)
  derived/
    model_input.csv                       — per-fish model input dataset (output of 02_prepare_model_input.ipynb)
    model_input_time_dependent.csv        — same but with time-dependent distances (output of S01)
    response_videos.csv                   — summary per response video
    all_positions.csv                     — 3D positions per fish per frame
    all_distances.csv                     — pairwise distances per video

notebooks/
  01_video_processing.ipynb              — processes raw stereo video annotations → 3D positions,
                                           distances, response categories
  02_prepare_model_input.ipynb           — filters observations and calculates neighbour influence inputs
  03_data_analysis.ipynb                 — statistical analyses (group size effects, spacing,
                                           initiator speed, speed decay parameter λ)
  04_simulations_analysis.ipynb          — analysis and comparison of simulation outputs vs empirical data
  05_figures.ipynb                       — all publication figures (Figures 1–4 + SI)
  S01_time_dependent_distances.ipynb     — supplementary: calculates time-dependent neighbour distances

scripts/
  04_model_fitting.R                     — fits all model variants (time decay × neighbour proximity ×
                                           speed × pooling), outputs AIC comparison table
  S01_model_fitting_time_dependent.R     — supplementary: fits models with time-dependent distances
  S02_model_fitting_different_D.R        — supplementary: fits models with alternative distance kernel
  simulations/
    00_config.R                          — simulation settings and paths (edit before running)
    01_io_and_helpers.R                  — I/O and helper functions
    02_load_model_fits.R                 — loads fitted model parameters
    03_fr_speed.R                        — FR speed distribution fitting
    04_empirical_loader.R                — loads empirical 3D networks
    04_gmm_loader.R                      — loads GMM-based 3D networks
    05_geometry.R                        — geometry and distance weight functions
    06_simulation_core.R                 — core simulation engine
    07_run_experiment.R                  — main simulation runner (entry point)
    08_merge_files.py                    — merges simulation output files
    09_trace_inputs.R                    — traces social inputs per fish over time
    10_calculate_Rt.R                    — calculates reproduction number R(t)

outputs/
  figures/                               — saved publication figures
  model_fitting/                         — model comparison tables, parameter estimates
  simulations/                           — simulation results (GMM_3D/, Real_3D/)
  data_analysis/                         — data analysis outputs (GMM values etc.)
```

---

## Model

The core model is a hazard-based social contagion model. The instantaneous response probability of fish *i* at time *t* is a function of:

- **Exponential time decay** since each neighbour responded (parameter λ)
- **Neighbour proximity** weighted by 1/*d*²
- **Initiator speed** raised to a fitted exponent *c*
- **Density term** capturing local neighbour proximity

Multiple model variants (combinations of time decay, proximity weighting, speed, and linear vs. threshold pooling) are compared by AIC in `04_model_fitting.R`. The best-performing model (DiTDeS with linear pooling) is used in all simulations.

---

## How to run

The pipeline runs in the following order. Processed data files are already provided in `data/`, so steps 1–2 can be skipped if you do not have access to the raw video annotations.

| Step | Script / Notebook | Description |
|------|-------------------|-------------|
| 1 | `notebooks/01_video_processing.ipynb` | Process stereo video annotations → 3D positions and response classifications |
| 2 | `notebooks/02_prepare_model_input.ipynb` | Filter observations and calculate neighbour influence inputs |
| 3 | `notebooks/03_data_analysis.ipynb` | Statistical analyses and estimation of speed decay parameter λ |
| 4 | `scripts/04_model_fitting.R` | Fit all model variants and produce AIC comparison table |
| 5 | `scripts/simulations/07_run_experiment.R` | Run simulations (edit `scripts/simulations/00_config.R` first to set paths and settings) |
| 6 | `scripts/simulations/09_trace_inputs.R` | Trace social inputs per fish over time |
| 7 | `scripts/simulations/10_calculate_Rt.R` | Calculate reproduction number R(t) |
| 8 | `notebooks/04_simulations_analysis.ipynb` | Analyse and compare simulation outputs against empirical data |
| 9 | `notebooks/05_figures.ipynb` | Generate all publication figures (Figures 1–4 + SI) |

Supplementary analyses (time-dependent distances, alternative distance kernel) are in `S01_*` and `S02_*` scripts and can be run independently after step 2.

---

## Dependencies

**Python** (notebooks):

- pandas, numpy, scipy, matplotlib, statsmodels, opencv-cv2, pathlib

**R** (scripts):

- dplyr, parallel, data.table, Matrix, here, tibble, readr

---

## Data note

Raw stereo video recordings are not included in this repository due to file size constraints. All processed data files required to reproduce the analyses, model fitting, simulations, and figures are provided in `data/` and `data/derived/`.

---

## Citation

This code will be linked to the published paper upon acceptance. Citation details will be added here at that time.
