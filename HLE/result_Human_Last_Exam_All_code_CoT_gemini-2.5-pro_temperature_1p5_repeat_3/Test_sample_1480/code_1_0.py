import numpy as np

# Given values for the firing levels
phi_k_x_prime = 0.7
mu_G_k_j_y_j = 0.9

# A t-norm operation is used to calculate the rule activation level.
# The most common t-norm is the minimum operator.
# Rule Activation Level = t_norm(phi_k_x_prime, mu_G_k_j_y_j)
rule_activation_level = np.minimum(phi_k_x_prime, mu_G_k_j_y_j)

# Print the final equation and the result
print(f"The rule activation level is calculated as the minimum (t-norm) of Phi_k(x') and μG_k_j(y_j).")
print(f"Rule Activation Level = min({phi_k_x_prime}, {mu_G_k_j_y_j})")
print(f"Rule Activation Level = {rule_activation_level}")
