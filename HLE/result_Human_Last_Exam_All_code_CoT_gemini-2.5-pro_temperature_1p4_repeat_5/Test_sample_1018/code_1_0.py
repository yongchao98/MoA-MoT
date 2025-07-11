import math

# Define the genera of the two surfaces
g1 = 31
g2 = 17

# The dimension of any surface is 2
n = 2
m = 2

# Calculate the simplicial volume of the first surface, ||Sigma_31||
# The formula for a surface with genus g > 1 is 4g - 4
sv1 = 4 * g1 - 4

# Calculate the simplicial volume of the second surface, ||Sigma_17||
sv2 = 4 * g2 - 4

# Calculate the binomial coefficient C(n+m, n)
# For this case, it is C(2+2, 2) = C(4, 2)
binomial_coeff = math.comb(n + m, n)

# Calculate the final simplicial volume of the product space
# ||Sigma_31 x Sigma_17|| = C(4, 2) * ||Sigma_31|| * ||Sigma_17||
total_sv = binomial_coeff * sv1 * sv2

# Print the final equation with all intermediate values
print(f"The simplicial volume of Sigma_31 x Sigma_17 is calculated using the product formula:")
print(f"||Sigma_31 x Sigma_17|| = C({n+m}, {n}) * ||Sigma_{g1}|| * ||Sigma_{g2}||")
print(f"First, ||Sigma_{g1}|| = 4*{g1} - 4 = {sv1}")
print(f"Second, ||Sigma_{g2}|| = 4*{g2} - 4 = {sv2}")
print(f"The binomial coefficient C({n+m}, {n}) is {binomial_coeff}.")
print(f"\nPutting it all together:")
print(f"Result = {binomial_coeff} * {sv1} * {sv2} = {total_sv}")