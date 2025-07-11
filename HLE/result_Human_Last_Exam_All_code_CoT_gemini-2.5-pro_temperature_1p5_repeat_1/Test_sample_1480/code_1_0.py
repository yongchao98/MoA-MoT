import numpy as np

# Step 1: Define the given values.
# Phi_k(x') is the firing level of the primary membership function.
phi_k_x_prime = 0.7

# μG_k_j(y_j) is the membership of the consequent parameter in the secondary fuzzy set.
mu_G_k_j = 0.9

# Step 2: Calculate the rule activation level using a t-norm.
# The most common t-norm is the minimum operator, T(a, b) = min(a, b).
# The activation level is the t-norm of the two given values.
activation_level = np.minimum(phi_k_x_prime, mu_G_k_j)

# Step 3: Print the explanation and the result.
print("The rule activation level for an interval type-3 fuzzy set is calculated by applying a t-norm to the primary firing level and the secondary membership value.")
print("We will use the minimum t-norm for this calculation.")
print(f"Given values: Phi_k(x') = {phi_k_x_prime}, μG_k_j(y_j) = {mu_G_k_j}")
print("\nEquation:")
print(f"Activation Level = min({phi_k_x_prime}, {mu_G_k_j})")
print(f"Activation Level = {activation_level}")