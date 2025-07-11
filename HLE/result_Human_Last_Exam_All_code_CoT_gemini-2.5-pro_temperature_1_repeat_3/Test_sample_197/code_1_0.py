import numpy as np

# Set a seed for reproducibility
np.random.seed(42)

# Number of data points in the observational study
n_samples = 10000

# Define the coefficients for the structural equations
# These represent the strength of the causal links.
# E -> A
c_ea = 0.8
# E -> D
c_ed = 0.7
# A -> B
c_ab = 0.6
# B -> C
c_bc = 0.5
# D -> C
c_dc = 0.4

# Generate data based on the structural equation model
# We assume the relationships are linear and add random noise.

# E is an exogenous variable (its cause is outside the model)
# We can simulate it from a standard normal distribution.
E = np.random.randn(n_samples)

# A is caused by E
noise_A = np.random.randn(n_samples)
A = c_ea * E + noise_A

# D is caused by E
noise_D = np.random.randn(n_samples)
D = c_ed * E + noise_D

# B is caused by A
noise_B = np.random.randn(n_samples)
B = c_ab * A + noise_B

# C is caused by B and D (C is a collider)
noise_C = np.random.randn(n_samples)
C = c_bc * B + c_dc * D + noise_C

# Now, let's test for correlation between A and D
# In our model, A does not cause D, and D does not cause A.
# However, they share a common cause E.
correlation_matrix = np.corrcoef(A, D)
correlation_ad = correlation_matrix[0, 1]

# The final equation is the calculation of the correlation coefficient.
# Let's print the numbers involved.
print(f"Number of samples: {n_samples}")
print(f"Coefficient for E->A (c_ea): {c_ea}")
print(f"Coefficient for E->D (c_ed): {c_ed}")
print(f"Calculated correlation between A and D: {correlation_ad:.4f}")
print("\nDespite no causal link between A and D, they are correlated because of the common cause E.")
print("Therefore, in this system, correlation does not imply causation.")
print("\nFinal Answer:")
print("No")
