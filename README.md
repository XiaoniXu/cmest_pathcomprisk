# Causal Mediation Analysis for Longitudinal Data with Competing Risks

## Overview
This repository contains the development, implementation, and validation of a comprehensive statistical pipeline for **Causal Mediation Analysis** in complex settings featuring **longitudinal mediators** and **competing risks**. 

It extends the state-of-the-art causal mediation framework by actively developing and validating the `cmcrhazard` R package—specifically the `cmest_pathcomprisk` command. This repository features robust methodologies to study complex aging-related health trajectories where patient outcomes track multi-timepoint physiological data and are constrained by dropouts or survival phenomena.

## Key Features
* **Time-Varying Mediation:** Accurately builds regression frameworks mapping continuous or time-varying variables across multiple measurement time points (e.g., baseline and follow-up mediators).
* **Survival & Competing Risk Integration:** Explicitly models time-to-event outcomes using Aalen's Additive Hazard Model (via the `timereg` package) combined with counting-process data structures.
* **Systematic Dropout Handling:** Uses sequential logistic regression models to robustly adjust for informative patient attrition or competing dropout events, preventing structural bias in long-term observational studies.
* **Exact-Match Bootstrapping (`refit = TRUE`):** Integrates parallel processing methodologies (`doParallel`) to execute computationally intensive, non-parametric resampling that yields highly precise standard errors and confidence intervals.

## Project Scope & Application
The methodologies and models developed in this repository were extensively validated against large-scale simulated distributions ($N = 7,500$ subjects, 1,000 bootstrap iterations) that mimic the multi-dimensional demographic and clinical profiles native to the **Strong Heart Study (SHS)**. 

### Core Tech Stack
* **Language:** R
* **Statistical Modeling:** `survival`, `timereg`, `mvtnorm`
* **Data Engineering:** `data.table`, `tidyr`, `dplyr`
* **Performance Optimization:** `doParallel`, `foreach`

## Selected Script Examples
* `test_med_longitudinal.R`: The core pipeline execution script demonstrating data staging, missingness tracking, survival modeling, and parallelized execution of the `cmest_pathcomprisk` module.
* `cmcrhazard/`: Contains the bespoke R package and utility logic actively driving the causal inference modeling framework.

---
### Author
**Justin Xiaoni Xu**
*Feel free to reach out via GitHub or connect with me via my personal website to discuss causal inference, biostatistics, or data science.*