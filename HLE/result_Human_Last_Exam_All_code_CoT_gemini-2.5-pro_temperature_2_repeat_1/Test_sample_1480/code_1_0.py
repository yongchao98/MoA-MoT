# The given values for the fuzzy rule
phi_k = 0.7  # Firing level of the antecedent: Phi_k(x′)
mu_G = 0.9   # Membership value of the consequent: μG_k_j(y_j)

# The t-norm operation combines these values. We will use the most common
# t-norm, which is the minimum operator.
# Rule Activation Level = T(Phi_k(x′), μG_k_j(y_j))
# Using minimum t-norm: Activation Level = min(Phi_k(x′), μG_k_j(y_j))

# Calculate the rule activation level
activation_level = min(phi_k, mu_G)

# Print the step-by-step calculation and the result
print(f"Calculating the rule activation level using the minimum t-norm:")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_G})")
print(f"The final rule activation level is: {activation_level}")
