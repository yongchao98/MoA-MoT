# The firing strength of the antecedent part of the fuzzy rule.
phi_k = 0.7

# The membership value of the consequent part of the fuzzy rule.
mu_g = 0.9

# A t-norm operation is used to combine these values to find the rule activation level.
# The most common t-norm is the minimum operation.
activation_level = min(phi_k, mu_g)

print("The rule activation level is calculated using the minimum t-norm operation:")
print(f"min({phi_k}, {mu_g}) = {activation_level}")