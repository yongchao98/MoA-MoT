import math

# Step 1: Define the genera and dimensions
g1 = 31
g2 = 17
dim_m = 2
dim_n = 2

# Step 2: Calculate the simplicial volume of the first surface, Sigma_31
# The formula for the simplicial volume of a surface of genus g >= 1 is ||Sigma_g|| = 4g - 4.
sv1 = 4 * g1 - 4

# Step 3: Calculate the simplicial volume of the second surface, Sigma_17
sv2 = 4 * g2 - 4

# Step 4: Calculate the binomial coefficient for the product formula
# The formula is C(m+n, m)
binom_coeff = math.comb(dim_m + dim_n, dim_m)

# Step 5: Calculate the final simplicial volume of the product space
total_sv = binom_coeff * sv1 * sv2

# Print the results, showing each number in the final equation
print(f"The simplicial volume of a surface Sigma_g (for g>=1) is given by the formula ||Sigma_g|| = 4g - 4.")
print(f"For the surface of genus g1 = {g1}, the simplicial volume is ||Sigma_{g1}|| = 4 * {g1} - 4 = {sv1}.")
print(f"For the surface of genus g2 = {g2}, the simplicial volume is ||Sigma_{g2}|| = 4 * {g2} - 4 = {sv2}.")
print(f"\nThe product formula for simplicial volumes is ||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||.")
print(f"For two surfaces (dimension 2), the binomial coefficient is C({dim_m}+{dim_n}, {dim_m}) = {binom_coeff}.")
print(f"\nTherefore, the simplicial volume of Sigma_{g1} x Sigma_{g2} is:")
print(f"||Sigma_{g1} x Sigma_{g2}|| = {binom_coeff} * {sv1} * {sv2} = {total_sv}")