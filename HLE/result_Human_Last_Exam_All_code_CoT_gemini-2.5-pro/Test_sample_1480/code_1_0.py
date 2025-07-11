# Given membership values for the antecedent parts of the fuzzy rule
phi_k = 0.7
mu_G_k_j = 0.9

# A t-norm is used to combine the membership values.
# The most common t-norm is the minimum operator.
activation_level = min(phi_k, mu_G_k_j)

# Print the equation and the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_G_k_j})")
print(f"Activation Level = {activation_level}")