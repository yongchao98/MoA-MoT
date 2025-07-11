import numpy as np

# Given values for the firing levels
phi_k = 0.7
mu_G_k_j = 0.9

# In fuzzy logic, the activation level of a rule is often determined by
# applying a t-norm operation to the antecedent and consequent firing levels.
# The most common t-norm is the minimum operator, T(a, b) = min(a, b).
# We will use the minimum t-norm to calculate the rule activation level.
activation_level = np.minimum(phi_k, mu_G_k_j)

# Print the calculation steps and the result
print("Calculating the rule activation level using the minimum t-norm:")
print(f"activation_level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"activation_level = min({phi_k}, {mu_G_k_j})")
print(f"The final rule activation level is: {activation_level}")