import matplotlib.pyplot as plt
import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

# --- Enzyme Kinetic Reaction Model ---
# Simplified denitrification reaction rate equation incorporating multiple substrate dependencies

def reaction_rate_denitrification(
    max_rate: float,  # μ_max_deni: Maximum specific reaction rate (1/s)
    biomass: float,   # X3: Microbial biomass concentration (g/L)
    biomass_saturation: float, # k_b3: Half-saturation constant for biomass (g/L)
    o2: float,        # O2: Dissolved oxygen concentration (mg/L)
    o2_inhibition: float, # k_O2I3: O2 inhibition constant (mg/L)
    n2o: float,       # N2O: Nitrous oxide concentration (mg/L)
    n2o_affinity: float, # K_N2O: N2O half-saturation constant (mg/L)
    formaldehyde: float, # CH2O: Formaldehyde concentration (mg/L)
    formaldehyde_affinity: float # K_CH2O: CH2O half-saturation constant (mg/L)
) -> float:
    """
    Calculates reaction rate considering:
    1. Biomass saturation kinetics
    2. Oxygen inhibition effect
    3. N2O consumption kinetics
    4. Formaldehyde consumption kinetics
    """
    biomass_term = biomass_saturation / (biomass_saturation + biomass)
    oxygen_inhibition = o2_inhibition / (o2_inhibition + o2)
    n2o_uptake = n2o / (n2o_affinity + n2o)
    formaldehyde_uptake = formaldehyde / (formaldehyde_affinity + formaldehyde)

    return max_rate * biomass * biomass_term * oxygen_inhibition * n2o_uptake * formaldehyde_uptake

# --- Experimental Simulation Setup ---
np.random.seed(42)
time_steps = 20  # 20 measurement time points

# Constant microbial biomass for simulation
BIOMASS = 0.3  # g/L (typical range depends on system)

# Simulated reactant profiles with measurement noise
O2_profile = np.linspace(3.0, 1.0, time_steps) + np.random.normal(0, 0.1, time_steps)  # mg/L
N2O_profile = np.linspace(1.2, 0.8, time_steps) + np.random.normal(0, 0.05, time_steps)  # mg/L
CH2O_profile = np.linspace(6.0, 3.0, time_steps) + np.random.normal(0, 0.2, time_steps)  # mg/L

# Ground truth kinetic parameters (hypothetical values for simulation)
TRUE_PARAMETERS = {
    "max_rate": 0.7,      # 1/s
    "biomass_saturation": 0.2,  # g/L
    "o2_inhibition": 0.3,  # mg/L
    "n2o_affinity": 0.4,   # mg/L
    "formaldehyde_affinity": 2.0  # mg/L
}

# Generate synthetic experimental data
true_rates = [
    reaction_rate_denitrification(
        TRUE_PARAMETERS["max_rate"],
        BIOMASS,
        TRUE_PARAMETERS["biomass_saturation"],
        O2_profile[t],
        TRUE_PARAMETERS["o2_inhibition"],
        N2O_profile[t],
        TRUE_PARAMETERS["n2o_affinity"],
        CH2O_profile[t],
        TRUE_PARAMETERS["formaldehyde_affinity"]
    )
    for t in range(time_steps)
]

observed_rates = np.clip(true_rates, a_min=0, a_max=None)  # Physical constraint: rates ≥ 0

# --- Bayesian Parameter Estimation Model ---
with pm.Model() as kinetic_model:

    # Priors based on biochemical knowledge
    max_rate = pm.Uniform('max_rate', lower=0.1, upper=1.5)  # Biologically plausible range
    biomass_saturation = pm.Uniform('biomass_saturation', lower=0.05, upper=1.0)
    o2_inhibition = pm.Uniform('o2_inhibition', lower=0.05, upper=1.0)
    n2o_affinity = pm.Uniform('n2o_affinity', lower=0.05, upper=1.0)
    formaldehyde_affinity = pm.Uniform('formaldehyde_affinity', lower=0.1, upper=5.0)

    # Measurement uncertainty models for reactants
    measured_o2 = pm.Normal('measured_o2', mu=O2_profile, sigma=0.05, shape=time_steps)
    measured_n2o = pm.Normal('measured_n2o', mu=N2O_profile, sigma=0.02, shape=time_steps)
    measured_ch2o = pm.Normal('measured_ch2o', mu=CH2O_profile, sigma=0.1, shape=time_steps)

    # Reaction rate prediction
    predicted_rate = pm.Deterministic(
        'predicted_rate',
        max_rate * BIOMASS *
        (biomass_saturation / (biomass_saturation + BIOMASS)) *
        (o2_inhibition / (o2_inhibition + measured_o2)) *
        (measured_n2o / (n2o_affinity + measured_n2o)) *
        (measured_ch2o / (formaldehyde_affinity + measured_ch2o))
    )

    # Measurement noise model
    rate_noise = pm.HalfNormal('rate_noise', sigma=0.01)  # Typical instrument precision
    rate_observations = pm.Normal(
        'rate_observations',
        mu=predicted_rate,
        sigma=rate_noise,
        observed=observed_rates
    )

    # MCMC sampling configuration
    inference_data = pm.sample(
        600,
        tune=600,
        target_accept=0.90,
        return_inferencedata=True
    )

# --- Results Analysis ---
# Parameter distributions
plt.figure(dpi=400)  # Increase DPI for sharper plot display
az.plot_posterior(
    inference_data,
    var_names=['max_rate', 'biomass_saturation', 'o2_inhibition',
               'n2o_affinity', 'formaldehyde_affinity'],
    figsize=(12, 8)
)
plt.suptitle('Posterior Distributions of Kinetic Parameters', y=1.02)
plt.tight_layout()
plt.show()
# plt.savefig("figures/posterior_distributions.png", dpi=300)  # Save with high DPI
# plt.show()


# Model validation
posterior_mean = inference_data.posterior['predicted_rate'].mean(dim=('chain', 'draw')).values
ci_lower = inference_data.posterior['predicted_rate'].quantile(0.05, dim=("chain", "draw"))
ci_upper = inference_data.posterior['predicted_rate'].quantile(0.95, dim=("chain", "draw"))

plt.figure(figsize=(10, 6), dpi=400)
plt.plot(observed_rates, 'o-', label='Observed Rates', lw=2)
plt.plot(posterior_mean, 's--', label='Model Prediction', lw=2)
plt.fill_between(range(time_steps), ci_lower, ci_upper, alpha=0.3, label='90% CI')
plt.xlabel('Time (arbitrary units)', fontsize=12)
plt.ylabel('Reaction Rate (units)', fontsize=12)
plt.title('Model Validation: Observed vs Predicted Reaction Rates', pad=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
