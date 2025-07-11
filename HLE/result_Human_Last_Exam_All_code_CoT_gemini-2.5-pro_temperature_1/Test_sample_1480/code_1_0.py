import numpy as np

# Given values for the rule activation
# Firing level of the antecedent part of the rule
phi_k = 0.7
# Membership value of the consequent part of the rule
mu_G_k_j = 0.9

# A t-norm is a type of binary operation used in fuzzy logic to model a fuzzy AND.
# The most common t-norm is the minimum operator, which we will use here.
activation_level = min(phi_k, mu_G_k_j)

print(f"The rule activation level is calculated by applying a t-norm (minimum operator) to the antecedent's firing level and the consequent's membership grade.")
print(f"Rule Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Rule Activation Level = min({phi_k}, {mu_G_k_j})")
print(f"The calculated rule activation level is: {activation_level}")