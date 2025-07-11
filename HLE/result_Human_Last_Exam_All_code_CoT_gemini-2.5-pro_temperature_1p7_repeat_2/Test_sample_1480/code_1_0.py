# Define the given membership values
phi_k_x_prime = 0.7
mu_G_k_j_y_j = 0.9

# The t-norm operation is a generalization of the logical AND operator.
# The most common t-norm is the minimum operator, which we will use here.
# Rule Activation Level = t_norm(Phi_k(x'), μG_k_j(y_j))
# Rule Activation Level = min(Phi_k(x'), μG_k_j(y_j))
activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

# Print the final equation and the result
print(f"The rule activation level using the minimum t-norm is calculated as:")
print(f"min({phi_k_x_prime}, {mu_G_k_j_y_j}) = {activation_level}")