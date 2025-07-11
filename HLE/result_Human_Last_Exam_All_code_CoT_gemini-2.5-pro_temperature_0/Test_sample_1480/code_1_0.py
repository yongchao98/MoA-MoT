# Define the given values
phi_k = 0.7
mu_G_k_j = 0.9

# The rule activation level is calculated using a t-norm operation.
# We will use the most common t-norm, which is the minimum operator.
activation_level = min(phi_k, mu_G_k_j)

# Print the final equation and the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_G_k_j}) = {activation_level}")