# Given values
phi_k_x_prime = 0.7
mu_G_k_j_y_j = 0.9

# The t-norm operation is typically the minimum function in this context.
# The rule activation level is calculated as the t-norm of the given values.
activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

# Print the calculation and the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"min({phi_k_x_prime}, {mu_G_k_j_y_j}) = {activation_level}")