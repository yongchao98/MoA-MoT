import numpy as np

# We found that the maximum ratio occurs for n=2.
# This script calculates the exact value of this ratio.

n = 2

# The radius r* which maximizes the classical distribution is r* = 3n^2 / 2
r_star = 3 * n**2 / 2.0

# Calculate the quantum radial distribution function D_n(r) for n=2 at r=r_star.
# D_2(r) = r^2 * [R_20(r)^2 + 3*R_21(r)^2]
# After simplification using the expressions for R_20 and R_21, we get:
# D_2(r) = (e^-r / 4) * (r^4 - 2r^3 + 2r^2)
Dn_r_star = (np.exp(-r_star) / 4.0) * (r_star**4 - 2 * r_star**3 + 2 * r_star**2)

# Calculate the classical radial distribution function D_n^c(r) at r=r_star.
# After normalization, this value is independent of n.
Dnc_r_star = 9 / (2 * np.pi * np.sqrt(3))

# Calculate the ratio
ratio = Dn_r_star / Dnc_r_star

# Output the numbers used in the final equation
print(f"The calculation for n={n}:")
print(f"Quantum distribution at r*={r_star}: D_n(r*) = {Dn_r_star}")
print(f"Classical distribution at r*={r_star}: D_n^c(r*) = {Dnc_r_star}")
print(f"Ratio = {Dn_r_star} / {Dnc_r_star}")

# Print the final result
print("The maximum value of the ratio is:")
print(ratio)