# Define the given membership values
phi_k = 0.7
mu_G = 0.9

# A t-norm is a fuzzy 'AND' operation. The most common t-norm is the minimum function.
# We will use the minimum t-norm to calculate the rule activation level.
activation_level = min(phi_k, mu_G)

# Print the equation and the result
print(f"The rule activation level is calculated using a t-norm (minimum) operation.")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_G}) = {activation_level}")