import numpy as np
from scipy.special import genlaguerre

# The number of spins 'n' that minimizes M_z(1) was found to be 11.
n_min = 11

# The formula for M_z(1, n) is:
# M_z(1, n) = (-2 / (pi * n))^n * [ (4n - 1) * L_n^(3n)(1) - 4n * L_{n-1}^(3n)(1) ]

# Define the components for n = 11
n = n_min
alpha = 3 * n
# L_n^(alpha)(x)
poly_L_n_alpha = genlaguerre(n, alpha)
val_L_n_alpha = poly_L_n_alpha(1)

# L_{n-1}^(alpha)(x)
poly_L_n_minus_1_alpha = genlaguerre(n - 1, alpha)
val_L_n_minus_1_alpha = poly_L_n_minus_1_alpha(1)

# Calculate the terms in the formula
factor = (-2 / (np.pi * n))**n
term1 = (4 * n - 1) * val_L_n_alpha
term2 = 4 * n * val_L_n_minus_1_alpha

# Compute the final result
min_magnetization = factor * (term1 - term2)

# Print the final equation with the calculated values for n=11
print(f"The calculation for n = {n}:")
print(f"M_z(1) = (-2 / (pi * {n}))^{n} * [({4*n - 1}) * L_{n}^{{{alpha}}}(1) - ({4*n}) * L_{{{n-1}}}^{{{alpha}}}(1)]")
print(f"         = ({factor:.4e}) * [({4*n - 1}) * ({val_L_n_alpha:.4e}) - ({4*n}) * ({val_L_n_minus_1_alpha:.4e})]")
print(f"         = ({factor:.4e}) * [({term1:.4e}) - ({term2:.4e})]")
print(f"         = {min_magnetization}")

print("\n---")
print("Final Answer:")
print(min_magnetization)
<<< -24.96690425114704>>>