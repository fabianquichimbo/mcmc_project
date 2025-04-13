#  Bayesian Inference of Denitrification Reaction Rate (r₃) using PyMC

This works presents a probabilistic modeling approach for estimating the denitrification reaction rate function \( r₃ \), using PyMC and time-series input data. The model incorporates uncertainty in kinetic parameters and key environmental drivers (O₂, N₂O, and CH₂O concentrations) to simulate and infer the reaction rate over time.

---

##  Reaction Function
The modeled denitrification rate function is defined as:
  
<img src="https://latex.codecogs.com/svg.image?\large&space;r_3=\mu_{\text{max}}^{\text{deni}}\cdot&space;X_3\cdot\left(\frac{k_{b3}}{k_{b3}&plus;X_3}\right)\cdot\left(\frac{k_{O2I3}}{k_{O2I3}&plus;O_2}\right)\cdot\left(\frac{N_2O}{K_{N_2O}&plus;N_2O}\right)\cdot\left(\frac{CH_2O}{K_{CH_2O}&plus;CH_2O}\right)" title="r_3=\mu_{\text{max}}^{\text{deni}}\cdot X_3\cdot\left(\frac{k_{b3}}{k_{b3}+X_3}\right)\cdot\left(\frac{k_{O2I3}}{k_{O2I3}+O_2}\right)\cdot\left(\frac{N_2O}{K_{N_2O}+N_2O}\right)\cdot\left(\frac{CH_2O}{K_{CH_2O}+CH_2O}\right)" />

Where:

  <img src="https://latex.codecogs.com/svg.image?\large&space;\(\mu_{\text{max}}^{\text{deni}}\)" title="\(\mu_{\text{max}}^{\text{deni}}\)" /> :  Maximum specific growth rate of denitrifiers. 

<img src="https://latex.codecogs.com/svg.image?\large&space;\;X_3\:" title="\;X_3\:" /> :Biomass concentration (assumed constant)

<img src="https://latex.codecogs.com/svg.image?\large&space;\;k_{b3},k_{O2I3},K_{N2O},K_{CH2O}\:" title="\;k_{b3},k_{O2I3},K_{N2O},K_{CH2O}\:" />: Kinetic constants

<img src="https://latex.codecogs.com/svg.image?\large&space;\(O_2,N_2O,CH_2O\):" title="\(O_2,N_2O,CH_2O\):" /> : Time-varying substrate concentrations

---

##  Data

Synthetic time-series data are generated with added Gaussian noise to emulate experimental conditions.
<p align="center">
  <img src="https://github.com/user-attachments/assets/50e45d3d-7b05-488f-9403-6ebf6579e5bb" alt="download (3)" />
</p>
<p align="center">
  <em>Figure 1: Rate production of Nitrou Oxide (N20) over time into a rock matrix.</em>
</p>

<p align="center">
  <em>Table 1: Representative theoretical kinetic parameters (Lee et al., 2011)</em>
</p>

| Parameter Name                | Symbol              | Value | Units |
|------------------------------|---------------------|-------|--------|
| Maximum denitrification rate | `μ_max`             | 0.7   | 1/s    |
| Biomass saturation constant  | `k_b3`              | 0.2   | g/L    |
| O₂ inhibition constant       | `k_O2I3`            | 0.3   | mg/L   |
| N₂O affinity constant        | `K_N2O`             | 0.4   | mg/L   |
| CH₂O (formaldehyde) affinity | `K_CH2O`            | 2.0   | mg/L   |

---

##  Model Description

- **Bayesian Framework** using [PyMC](https://www.pymc.io/)

This code section implements a Bayesian parameter estimation model for biochemical reaction kinetics using PyMC. It defines uniform prior distributions for five kinetic parameters (max_rate, biomass_saturation, o2_inhibition, n2o_affinity, formaldehyde_affinity) based on biologically plausible ranges1, models measurement uncertainty for one reactant (N₂O) via normal distributions centered around experimental profiles, and constructs a deterministic reaction rate equation incorporating saturation/inhibition terms. The model accounts for observation noise through a half-normal distributed rate_noise term and uses Monte Carlo sampling (600 tuning + 600 sampling steps, 90% target acceptance) to infer posterior distributions. The pm.Deterministic node explicitly tracks the predicted reaction rate for analysis, while observed data (observed_rates) constrains the parameter estimation through the likelihood function.

Note: Uniform priors reflect initial uncertainty about parameters while respecting biochemical constraints. The rate equation structure follows common enzyme kinetics formulations used in systems biology. MCMC configuration balances computational cost with sampling quality for parameter posteriors, and Bayesian updating combines prior knowledge with experimental data through the likelihood.

---

### Posterior Distributions

<p align="center">
  <img src="https://github.com/user-attachments/assets/1b7f58d5-c79a-41f6-8a42-5b62ebcfc7b0" alt="download (1)" />
</p>
<p align="center">
  <em>Figure 2: posterior distributions of the kinetic parameters estimated using Markov Chain Monte Carlo (MCMC) sampling.</em>
</p>

In general, this figure displays the posterior distributions of the kinetic parameters estimated using Markov Chain Monte Carlo (MCMC) sampling. Notably, the "n2o_affinity" plot shows the distribution of the affinity constant for nitrous oxide (N₂O) in the modeled biological system, with a mean value of approximately 0.48. Each plot illustrates the uncertainty in the estimated parameter value after incorporating the observed data, with the 94% highest density interval (HDI) providing a range for each parameter. The shape and spread of the N₂O affinity distribution, along with those of other parameters, reflect the sensitivity of the model to these parameters and the information content of the data. The HDI represents the range within which the parameter value is most likely to fall, given the model and the data and the narrower distributions indicate more precise parameter estimates, while broader distributions suggest greater uncertainty.

### Predicted vs Observed \( r_3 \)

<p align="center">
  <img src="https://github.com/user-attachments/assets/4ba15f13-e925-49e2-843f-e4035fe27a45" alt="download (2)" />
</p>
<p align="center">
  <em>Figure 3: Assessing of how well the Bayesian model informed by MCMC sampling can captures the current dynamics of the system.</em>
</p>

On figure 3, on the x-axis of "Time" indicates the progression of the reaction or process being modeled. The y-axis, labeled "Reaction Rate" quantifies the speed at which the reaction proceeds. The plot overlays two key datasets: "Observed Rates" (blue line with circular markers): These are the "experimentally" measured reaction rates at different time points. They represent the real-world data that the model aims to replicate. "Model Prediction" (orange line with square markers): This line represents the reaction rates predicted by the Bayesian model, using the kinetic parameters estimated via MCMC. 

The light blue shaded area around the "Model Prediction" line represents the 90% Credible Interval (CI). This interval signifies the range within which the model predicts the true reaction rate to lie with 90% probability, given the estimated parameter distributions2. A narrow credible interval indicates higher confidence in the model's predictions. By visually comparing the "Observed Rates" with the "Model Prediction" and considering the width of the 90% CI, one can assess the goodness-of-fit of the model. If the "Observed Rates" closely follow the "Model Prediction" and fall within the 90% CI, it suggests that the model accurately captures the underlying process. Conversely, significant deviations between the observed data and the model's predictions, or a wide credible interval, may indicate model misspecification or the need for additional data.

Model validation is an essential step to ensure that the model is an accurate representation of the real system. The credible interval reflects the uncertainty in the model predictions due to the uncertainty in the estimated parameters, and a good fit between the observed and predicted rates provides confidence in the model's ability to make accurate predictions.

---



