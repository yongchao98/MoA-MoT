# Define the given values
phi_k = 0.7
mu_g = 0.9

# The rule activation level (w_k) is calculated using a t-norm.
# The most common t-norm is the minimum operator.
activation_level = min(phi_k, mu_g)

# Print the explanation and the final equation with the result
print("The rule activation level w_k(x') is calculated using the minimum t-norm as follows:")
print(f"w_k(x') = min(Phi_k(x'), μG_k_j(y_j))")
print(f"w_k(x') = min({phi_k}, {mu_g}) = {activation_level}")