import numpy as np

# --- Given Values ---
# Phi_k(x′) represents the firing level of the antecedent part of the rule.
phi = 0.7
# muG_k_j(y_j) represents a membership value from the consequent part of the rule.
mu_G = 0.9

# --- T-norm Operation ---
# A t-norm is a function that models the logical AND operation in fuzzy logic.
# The most common t-norm is the minimum function, t(a, b) = min(a, b).
# We calculate the rule activation level by applying this t-norm to the given values.
activation_level = min(phi, mu_G)

# --- Output the Result ---
print("Calculating the rule activation level using the minimum t-norm.")
print(f"The given values are Phi_k(x') = {phi} and mu_G_k_j(y_j) = {mu_G}.")
print(f"The calculation is: min({phi}, {mu_G})")
print(f"The resulting rule activation level is: {activation_level}")
