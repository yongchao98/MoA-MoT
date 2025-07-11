import numpy as np

# Given values
phi_k = 0.7  # Firing level of the antecedent part
mu_G_k_j = 0.9  # Membership degree of the consequent part

# The rule activation level is computed by applying a t-norm to the antecedent
# firing level and the consequent membership function.
# We will use the minimum t-norm, which is the most common.
# t_norm(a, b) = min(a, b)
activation_level = np.minimum(phi_k, mu_G_k_j)

print("The rule activation level is the result of the t-norm operation between the antecedent firing level and the consequent's membership degree.")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_G_k_j})")
print(f"The calculated rule activation level is: {activation_level}")