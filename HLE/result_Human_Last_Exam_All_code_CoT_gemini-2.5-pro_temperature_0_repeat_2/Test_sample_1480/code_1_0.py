# Define the given values
# Phi_k(x') is the firing level of the rule's antecedent
phi_k = 0.7
# muG_k_j(y_j) is the membership degree of the consequent
mu_g = 0.9

# The rule activation level is calculated by applying a t-norm operation
# to the antecedent's firing level and the consequent's membership degree.
# We will use the minimum t-norm, which is the most common choice.
activation_level = min(phi_k, mu_g)

# Print the final equation showing the inputs and the result
print(f"The rule activation level is calculated as: min({phi_k}, {mu_g})")
print(f"Result: {activation_level}")