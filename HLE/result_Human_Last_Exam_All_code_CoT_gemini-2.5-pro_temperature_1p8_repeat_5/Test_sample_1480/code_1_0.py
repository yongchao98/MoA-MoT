# Define the given values for the fuzzy set parameters
phi_k_x_prime = 0.7
mu_G_k_j_y_j = 0.9

# A t-norm is a fuzzy logic operator for intersection (AND).
# The most commonly used t-norm is the minimum function.
# We calculate the rule activation level by taking the minimum
# of the antecedent firing level and the consequent membership grade.
activation_level = min(phi_k_x_prime, mu_G_k_j_y_j)

# Print the final equation with each number and the result
print(f"The rule activation level is calculated using the minimum t-norm operation.")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k_x_prime}, {mu_G_k_j_y_j})")
print(f"Result = {activation_level}")