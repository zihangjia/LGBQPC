# LGBQPC: Local Granular-Ball Quality Peaks Clustering

This repository contains the code and datasets used in the paper, which is submitted to IEEE TCYB:

**LGBQPC: Local Granular-Ball Quality Peaks Clustering**  

---

## Overview

The **LGBQPC** algorithm is designed for clustering tasks, particularly on datasets with complex structures or non-uniform densities. This repository includes all code and small-scale datasets used in our experiments on D1–D30 datasets.

---

## Contents

- `LGBQPC/` : MATLAB implementation of the LGBQPC algorithm.
- `datasets/` : All small-scale datasets (D1–D30) used in the experiments.
- `README.md` : This file.

---

## Requirements

- **MATLAB 2024b** : Required to run the LGBQPC algorithm.
- **Python** : Required for evaluating clustering results (metrics computation, etc.).

---

## Usage

1. Open MATLAB and navigate to the `LGBQPC` folder.
2. Run the `main.m` function to execute the LGBQPC algorithm on datasets D1–D30:

```matlab
main
