import numpy as np

# Given membership values
phi_k = 0.7
mu_G_k_j = 0.9

# A t-norm operation is a fuzzy logic AND. The most common t-norm is the minimum operator.
# We will calculate the rule activation level using the minimum t-norm.
activation_level = np.minimum(phi_k, mu_G_k_j)

# Print the final equation with all the numbers.
print(f"The rule activation level (τ) is calculated using the minimum t-norm:")
print(f"τ = min(Phi_k(x'), μG_k_j(y_j))")
print(f"τ = min({phi_k}, {mu_G_k_j}) = {activation_level}")