# Code for "Guidance on Individualized Treatment Rule Estimation in High Dimensions"

> [Philippe Boileau](https://pboileau.ca/), Ning Leng and [Sandrine
> Dudoit](https://www.stat.berkeley.edu/~sandrine/)

This repository contains the code required to reproduce, and even extend, the
simulation study presented in the paper. These simulations rely on the
[`simChef`](https://github.com/Yu-Group/simChef) `R` package's flexible
framework.

Here's a brief description of the repository's contents:

- `R/`: Where all of the `R` code is saved and organized in the following
  sub-directories.
  - `R/dgp-functions/`: Data-generating process functions.
  - `R/method-functions/`: Functions implementing (or wrapping existing
    implementations of) the considered methods.
  - `R/eval-functions/`: Functions for summarizing the simulation study results.
  - `R/viz-functions/`: Result visualization functions.
  - `R/meals/`: Scripts defining the various simulation experiments.
- `hpc/`: An empty folder for organizing scripts that run the contents of
  `R/meals/` in a high performance computing environment (HPC). Doing so is
  recommended over reproducing these simulations locally. Bash files for HPCs
  using the SLURM workload manager are available upon request.
- `logs/`: An empty folder for collecting messages output by the HPC.
- `results`: An empty folder for collecting simulation study results.
