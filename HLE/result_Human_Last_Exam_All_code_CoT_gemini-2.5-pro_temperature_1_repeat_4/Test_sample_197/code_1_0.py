import numpy as np
from scipy.stats import pearsonr

# Set a seed for reproducibility of the simulation
np.random.seed(42)

# Number of data points to simulate
n_samples = 10000

# Define the causal coefficients for the structural equations
# These are the numbers that define the strength of the causal links
coef_E_A = 0.8  # Strength of E -> A
coef_E_D = 0.7  # Strength of E -> D
coef_A_B = 0.6  # Strength of A -> B
coef_B_C = 0.5  # Strength of B -> C
coef_D_C = 0.4  # Strength of D -> C

# --- Data Simulation ---
# We generate data based on the structural equations of the model.
# E is the common cause, which we model as a random variable.
E = np.random.randn(n_samples)

# A is caused by E, following the equation: A = 0.8 * E + noise
A = coef_E_A * E + np.random.randn(n_samples)

# D is also caused by E, following the equation: D = 0.7 * E + noise
D = coef_E_D * E + np.random.randn(n_samples)

# B is caused by A (not needed for the A-D correlation, but included for model completeness)
B = coef_A_B * A + np.random.randn(n_samples)

# C is caused by B and D (C is a "collider")
C = coef_B_C * B + coef_D_C * D + np.random.randn(n_samples)


# --- Analysis ---
# We now test for correlation between A and D. In our model, A does not cause D,
# and D does not cause A. The question is whether the correlation we might find
# implies a causal relationship.

# Calculate the Pearson correlation coefficient and the p-value.
corr_AD, p_value_AD = pearsonr(A, D)

print("Analysis of the Structural Equation Model: E->A->B->C<-D<-E")
print("\nWe are testing if an observed correlation between A and D implies causation.")
print("\nStructural Equations used for simulation:")
print(f"A = {coef_E_A} * E + noise")
print(f"D = {coef_E_D} * E + noise")
print("\nNote: There is no direct causal term between A and D in these equations.")

print("\n--- Simulation Result ---")
# The final equation is the calculation of the correlation coefficient. We output its result.
print(f"The calculated Pearson correlation between A and D is: {corr_AD:.4f}")
print(f"The p-value is: {p_value_AD:.2e}")

print("\n--- Conclusion ---")
print("The simulation shows a statistically significant correlation between A and D.")
print("This correlation exists because both A and D are caused by the common variable E (a confounder).")
print("It does not exist because of a causal link between A and D.")
print("\nTherefore, in this system, correlation does not imply causation.")
