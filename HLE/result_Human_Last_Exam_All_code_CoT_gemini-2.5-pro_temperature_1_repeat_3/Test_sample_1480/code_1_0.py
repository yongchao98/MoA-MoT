import numpy as np

# Given values for the fuzzy set
phi_k = 0.7
mu_G = 0.9

# A t-norm is a type of binary operation used in fuzzy logic to model a logical conjunction (AND).
# The most common t-norm is the minimum operator, which we will use here.
# The rule activation level is calculated by applying the t-norm to the given values.
activation_level = np.minimum(phi_k, mu_G)

# Print the calculation and the result
print(f"The rule activation level is calculated using the minimum t-norm:")
print(f"min({phi_k}, {mu_G}) = {activation_level}")