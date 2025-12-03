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
