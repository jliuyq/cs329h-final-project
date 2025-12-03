# cs329h-final-project
## Overview
This repository contains the R implementation for the final project simulation study. It implements an Adaptive Pairwise Comparison framework using Bradley-Terry models. The goal is to evaluate different active learning policies (Random, Entropy, and CAT) on their ability to recover item parameters under different clustering conditions.

## Repository Structure
* `main.R`: The primary script containing all simulation logic, model definitions, and experiment loops.
* `README.md`: Project documentation.
* `bt_cat_summary_metrics.csv`: (Generated Output) Contains the aggregated metrics from the simulation.

## Requirements & Setup
The project is built using R (version 4.0+ recommended).

### Dependencies
The simulation relies on base R and the `dplyr` package for data aggregation.

To install dependencies, run:
```R
install.packages("dplyr")
```

## User Guide: How to Reproduce Results
### Step 1: File Setup

Ensure that `main.R` is in your current working directory. You can check your working directory in R using `getwd()` and set it using `setwd("/path/to/your/folder")`.

### Step 2: Full Experiment (Replicate Paper) To reproduce the full results presented in the project report:
* Open `main.R`.
* Scroll to Section 5 (Main Execution Block) at the bottom of the script.
* Modify the `exp_config` to the values you like.
* Run the script. Note: This may take 15–30 minutes depending on your hardware.
