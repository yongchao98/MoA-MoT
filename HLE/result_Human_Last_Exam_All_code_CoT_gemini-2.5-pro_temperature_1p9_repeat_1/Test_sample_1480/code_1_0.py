import numpy as np

# Step 1: Define the given values.
# Phi_k(x') is the firing strength of the antecedent part of the rule.
phi_k = 0.7
# muG_k_j(y_j) is the membership degree in the consequent fuzzy set.
mu_g = 0.9

# Step 2: Apply the t-norm operation.
# We will use the minimum t-norm, as it's the most common for fuzzy implication.
# The activation level is the result of clipping the consequent membership by the antecedent firing strength.
activation_level = np.minimum(phi_k, mu_g)

# Step 3: Print the result in a descriptive equation.
# The final output shows each number used in the calculation.
print(f"The rule activation level is found by applying a t-norm (minimum) to the antecedent firing strength and the consequent membership degree.")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_g}) = {activation_level}")