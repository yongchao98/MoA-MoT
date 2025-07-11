# Given values
phi_k = 0.7
mu_G_k_j = 0.9

# A t-norm operation is a generalized intersection.
# The most common t-norm is the minimum function.
# We will use the minimum t-norm to calculate the rule activation level.
activation_level = min(phi_k, mu_G_k_j)

# Print the final equation and the result
print(f"The rule activation level is calculated using the minimum t-norm operation.")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_G_k_j}) = {activation_level}")