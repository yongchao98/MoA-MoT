# Given values for the fuzzy rule
phi_k = 0.7  # Firing level of the antecedent
mu_G_k_j = 0.9 # Membership grade of the consequent

# The rule activation level is calculated using a t-norm operation.
# We will use the most common t-norm, which is the minimum operator.
activation_level = min(phi_k, mu_G_k_j)

# Print the calculation steps and the final result
print("The rule activation level is the t-norm of the antecedent's firing level and the consequent's membership grade.")
print(f"Using the minimum t-norm, the calculation is:")
print(f"Rule Activation Level = min({phi_k}, {mu_G_k_j})")
print(f"Rule Activation Level = {activation_level}")