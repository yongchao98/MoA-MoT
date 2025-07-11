# Define the given values for the rule
# Firing level of the antecedent part of the rule
phi_k = 0.7
# Membership grade of the consequent part of the rule
mu_G = 0.9

# The rule activation level is determined by applying a t-norm operation
# to the antecedent's firing level and the consequent's membership grade.
# The most common t-norm is the minimum (min) operator.
# Activation Level = t-norm(phi_k, mu_G)
activation_level = min(phi_k, mu_G)

# Print the calculation steps and the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"Activation Level = min({phi_k}, {mu_G}) = {activation_level}")