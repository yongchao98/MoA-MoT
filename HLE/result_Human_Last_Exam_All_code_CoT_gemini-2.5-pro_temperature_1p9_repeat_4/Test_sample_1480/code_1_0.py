import numpy as np

# Step 1: Define the antecedent firing level.
phi_k = 0.7

# Step 2: Define the consequent membership grade.
mu_g = 0.9

# Step 3: Apply the t-norm operation.
# The most common t-norm is the minimum operator, which we use here.
# This calculates the activation of the rule's consequent.
activation_level = np.minimum(phi_k, mu_g)

# Step 4: Print the final equation and the result.
print(f"The rule activation level is calculated as the t-norm (minimum) of the antecedent firing level and the consequent membership grade.")
print(f"Activation Level = min(Phi_k(x'), μG_k_j(y_j))")
print(f"Activation Level = min({phi_k}, {mu_g}) = {activation_level}")