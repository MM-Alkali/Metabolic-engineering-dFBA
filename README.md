# dFBA: Dynamic Flux Balance Analysis of Pyruvate Routing in Yeast

Dynamic Flux Balance Analysis of pyruvate routing in *Saccharomyces cerevisiae* under aerobic and anaerobic conditions with CDC19 (pyruvate kinase) constraints. MATLAB pipeline for simulating short-term metabolic responses to 150‑second glucose pulses, with comprehensive sensitivity and flux variability analyses.

## Overview

This repository contains the complete code and data for dynamic Flux Balance Analysis (dFBA) simulations of *Saccharomyces cerevisiae* central carbon metabolism under glucose pulse conditions. The study investigates **oxygen-dependent pyruvate routing** under CDC19 (pyruvate kinase) constraints, comparing single-pulse (SP) and repetitive-pulse (RP) conditions under both aerobic and anaerobic environments.

**Key features:**
- MATLAB pipeline for time-resolved dFBA simulations
- Four experimental conditions: Aerobic_SP, Aerobic_RP, Anaerobic_SP, Anaerobic_RP
- Pyruvate kinase (CDC19) constraints: **60% of theoretical capacity under anaerobic** (feed phase), **unconstrained under aerobic**
- 8 main metabolic parameters including glucose uptake, growth rate, CO₂, O₂, ethanol, pyruvate, PYK, and PDC/PDH routing
- Pyruvate routing analysis (3x4 layout) showing PYK, PDC, and PDH flux profiles
- Systematic sensitivity analysis (0–100% PYK capacity) with automatic optimal constraint identification
- Flux Variability Analysis (FVA) to assess metabolic flexibility at 95% optimal growth
- Model validation figures comparing predictions to experimental data
- Publication‑ready figures with separate legend panels

## Results Summary

| Condition | PYK Constraint | PYK Flux (mmol/gDW/h) | PDC/PYK Ratio | Ethanol Yield | Growth Rate |
|-----------|---------------|----------------------|---------------|---------------|-------------|
| Anaerobic_SP | 60% | 57–75 | 100% | 2.0 | ∼0 |
| Anaerobic_RP | 60% | 64–66 | 100% | 2.0 | ∼0 |
| Aerobic_SP | Unconstrained | ∼13 | 5% | 0.0 | 0.45 |
| Aerobic_RP | Unconstrained | ∼13 | 5% | 0.0 | 0.45 |


## 🔧 Requirements

### MATLAB (R2023b or later)
- [COBRA Toolbox v3.0](https://opencobra.github.io/)
- [GLPK](https://www.gnu.org/software/glpk/) (or any supported LP solver)
- Statistics and Machine Learning Toolbox (for table operations)

## 🚀 Usage

1. Clone the repository
2. Place experimental data Excel files in the data directory
3. Run `Run.m` in MATLAB

The script will:
- Load the yeast-GEM model
- Run dFBA for all four conditions
- Generate main comparison figures (3x4 layout)
- Perform sensitivity analysis (0–100% PYK capacity)
- Perform Flux Variability Analysis (FVA)
- Save all figures and CSV results to the output directory

## 📊 Output Figures

| Figure | Description |
|--------|-------------|
| `Main_Comparison_SP_RP.png` | 12‑panel figure (8 main parameters + 4 pyruvate routing panels) |
| `PYK_Sensitivity_Aerobic.png` | Sensitivity analysis for aerobic conditions (2x4 layout) |
| `PYK_Sensitivity_Anaerobic.png` | Sensitivity analysis for anaerobic conditions (2x4 layout) |
| `FVA_Combined.png` | Flux Variability Analysis results |
| `Supplementary_Model_Validation_*.png` | Model vs experimental validation per condition |

## 📝 Citation

If you use this code, please cite:

Muawiya, M.A et al. (2026)...

## 📧 Contact

[muawiyamuhd@mail.ecust.edu.cn]
