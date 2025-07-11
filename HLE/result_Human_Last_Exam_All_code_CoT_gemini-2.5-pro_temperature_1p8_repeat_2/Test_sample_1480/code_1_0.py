# Given values for the fuzzy set rule
phi_k_x_prime = 0.7
mu_G_k_j_y_j = 0.9

# The t-norm operation is commonly implemented using the minimum function
# to calculate the rule activation level (firing strength).
activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

# Print the calculation and the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"min({phi_k_x_prime}, {mu_G_k_j_y_j}) = {activation_level}")