import numpy as np

# Given membership values
phi_k = 0.7
mu_g = 0.9

# The t-norm operation combines the membership values. We will use the most
# common t-norm, which is the minimum operator.
activation_level = np.minimum(phi_k, mu_g)

# Print the final equation and the result
print("The rule activation level (w) is calculated using the minimum t-norm operation:")
print(f"w = min({phi_k}, {mu_g}) = {activation_level}")