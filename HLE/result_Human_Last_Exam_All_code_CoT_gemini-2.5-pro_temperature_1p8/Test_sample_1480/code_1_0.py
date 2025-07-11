import numpy as np

# Define the given values for the k-th rule
# Firing level of the antecedent part
phi_k = 0.7
# Membership of the j-th consequent
mu_g = 0.9

# The rule activation level is calculated using a t-norm operation.
# We will use the minimum t-norm, which is the most common choice.
# Activation = t-norm(phi_k, mu_g)
rule_activation = np.minimum(phi_k, mu_g)

# Print the result and the equation
print("The rule activation level is the result of the t-norm operation on the antecedent's firing level and the consequent's membership degree.")
print("Using the minimum t-norm:")
print(f"Activation = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation = min({phi_k}, {mu_g})")
print(f"The calculated rule activation level is: {rule_activation}")