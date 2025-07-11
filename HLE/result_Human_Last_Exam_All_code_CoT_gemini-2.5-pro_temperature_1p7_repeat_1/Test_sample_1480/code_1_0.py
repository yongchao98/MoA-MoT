# Define the given membership values
phi_k_x_prime = 0.7
mu_G_k_j_y_j = 0.9

# A t-norm operation is a form of fuzzy 'AND'. The most common t-norm is the
# minimum operator, which we will use here to find the rule activation level.
rule_activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

# Print the final equation with the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"min({phi_k_x_prime}, {mu_G_k_j_y_j}) = {rule_activation_level}")