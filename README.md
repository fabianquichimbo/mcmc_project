#  Bayesian Inference of Denitrification Reaction Rate (r₃) using PyMC

This repository presents a probabilistic modeling approach for estimating the denitrification reaction rate function \( r_3 \), using PyMC and time-series input data. The model incorporates uncertainty in kinetic parameters and key environmental drivers (O₂, N₂O, and CH₂O concentrations) to simulate and infer the reaction rate over time.

---

##  Reaction Function

The modeled denitrification rate function is defined as:

\[
r_3 = \mu_{\text{max}}^{\text{deni}} \cdot X_3 \cdot \left(\frac{k_{b3}}{k_{b3} + X_3}\right) \cdot \left(\frac{k_{O2I3}}{k_{O2I3} + O_2}\right) \cdot \left(\frac{N_2O}{K_{N_2O} + N_2O}\right) \cdot \left(\frac{CH_2O}{K_{CH_2O} + CH_2O}\right)
\]

Where:
- \( \mu_{\text{max}}^{\text{deni}} \): Maximum specific growth rate of denitrifiers  
- \( X_3 \): Biomass concentration (assumed constant)  
- \( k_{b3}, k_{O2I3}, K_{N2O}, K_{CH2O} \): Kinetic constants  
- \( O_2, N_2O, CH_2O \): Time-varying substrate concentrations

---

##  Simulated Data

Synthetic time-series data are generated with added Gaussian noise to emulate experimental conditions.

<div align="center">
  <img src="figures/input_data.png" alt="Input Time-Series Data" width="600"/>
</div>

---

##  Model Description

- **Bayesian Framework** using [PyMC](https://www.pymc.io/)
- Uncertain inputs modeled with normal priors (O₂, N₂O, CH₂O)
- Kinetic parameters modeled with uniform priors
- Observed \( r_3 \) includes noise modeled as a normal likelihood with HalfNormal prior on standard deviation

---

##  Inference Results

### Posterior Distributions

<div align="center">
  <img src="figures/posterior_distributions.png" alt="Posterior Distributions" width="700"/>
</div>

### Predicted vs Observed \( r_3 \)

<div align="center">
  <img src="figures/predicted_vs_observed.png" alt="Predicted vs Observed r3" width="700"/>
</div>

---

##  Project Structure


