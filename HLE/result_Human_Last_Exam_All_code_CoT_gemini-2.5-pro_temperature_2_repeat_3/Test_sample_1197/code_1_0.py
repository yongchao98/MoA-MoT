import numpy as np

# Set a seed for consistent, reproducible results.
np.random.seed(42)

print("Modeling Challenge to Genetic Drift Predictions")
print("="*50)
print("The predictive models of genetic drift often use the synonymous substitution rate (dS) as a neutral baseline.")
print("Deviations in the nonsynonymous rate (dN) from this baseline are then attributed to selection.")
print("The challenge arises when dS and dN are not independent.")

print("\n--- Scenario 1: Ideal Model (Independent Rates) ---")
print("In this model, dN and dS are independent. A high dN/dS ratio would clearly suggest positive selection.")

# Simulate 3 genes. dS is a proxy for the mutation/drift rate.
ds_independent = np.random.normal(loc=0.020, scale=0.003, size=3)
# dN is very low due to purifying selection and independent of dS.
dn_independent = np.random.normal(loc=0.004, scale=0.001, size=3)

# Ensure no negative values for rates
ds_independent[ds_independent < 0] = 0.001
dn_independent[dn_independent < 0] = 0.001

for i in range(3):
    rate_equation = dn_independent[i] / ds_independent[i]
    # Here we print each number in the final equation as requested
    print(f"Gene {i+1} Equation: dN({dn_independent[i]:.5f}) / dS({ds_independent[i]:.5f}) = {rate_equation:.3f}")

print("\n--- Scenario 2: Correlated Rates (The Challenge) ---")
print("Here, regional mutation biases cause dS and dN to be correlated. This complicates interpretation.")
print("A high observed dN might just reflect a high local mutation rate, not positive selection.")

# Simulate 3 genes where dN is correlated with dS.
ds_correlated = np.random.normal(loc=0.020, scale=0.008, size=3) # Use a larger scale for more variation
# dN is now a function of dS (e.g., 0.2 * dS) plus some random noise.
correlation_factor = 0.2
noise = np.random.normal(loc=0, scale=0.0005, size=3)
dn_correlated = (ds_correlated * correlation_factor) + noise

# Ensure no negative values for rates
ds_correlated[ds_correlated < 0] = 0.001
dn_correlated[dn_correlated < 0] = 0.001

for i in range(3):
    rate_equation = dn_correlated[i] / ds_correlated[i]
    # For Gene 3, both dN and dS are high, but the ratio is not elevated.
    # Without knowing about the correlation, one might misinterpret the high dN value.
    # This shows the challenge described in Choice B.
    print(f"Gene {i+1+3} Equation: dN({dn_correlated[i]:.5f}) / dS({ds_correlated[i]:.5f}) = {rate_equation:.3f}")

print("\nConclusion: The simulation demonstrates that a correlation between dS and dN makes dS an unreliable")
print("baseline for drift, which directly challenges the predictive models.")