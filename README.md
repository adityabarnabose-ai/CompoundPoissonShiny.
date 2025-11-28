# Compound Poisson Process Simulator

This repository contains an **R Shiny app** to simulate a **Compound Poisson Process** with exponential inter-arrival times and exponential jump sizes. It also includes a **LaTeX report** with theory, derivations, and plots.

---

## Features

- Simulates the **number of arrivals N(t)** with a Poisson process.
- Simulates **jump sizes X_i** using an exponential distribution.
- Calculates and plots the **cumulative sum S(t)** (the compound process).
- Generates a **histogram of jump sizes**.
- Interactive inputs for:
  - Arrival rate (λ)
  - Jump rate (μ)
  - Simulation horizon (t)

---

## Requirements

- **R version 4.5 or higher**  
- R packages:
  - `shiny`
  - `ggplot2`

Install packages using:

```r
install.packages("shiny")
install.packages("ggplot2")
