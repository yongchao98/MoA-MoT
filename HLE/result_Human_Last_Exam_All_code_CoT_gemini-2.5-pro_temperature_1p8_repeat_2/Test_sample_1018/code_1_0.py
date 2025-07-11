import math

# Define the genera of the surfaces
g1 = 31
g2 = 17

# The dimension of a surface is 2
m = 2
n = 2

# Step 1: Calculate the simplicial volume of Sigma_g1
# The formula for g >= 2 is ||Sigma_g|| = 4 * (g - 1)
sv1 = 4 * (g1 - 1)

# Step 2: Calculate the simplicial volume of Sigma_g2
sv2 = 4 * (g2 - 1)

# Step 3: Calculate the binomial coefficient for the product formula
# C(m+n, m)
total_dim = m + n
coeff = math.comb(total_dim, m)

# Step 4: Calculate the final simplicial volume of the product
final_sv = coeff * sv1 * sv2

# Print the result in the format of the full equation
print(f"The simplicial volume of Σ_{g1} x Σ_{g2} is computed using the product formula:")
print(f"||Σ_{g1} x Σ_{g2}|| = C({total_dim}, {m}) * ||Σ_{g1}|| * ||Σ_{g2}||")
print(f"||Σ_{g1}|| = 4 * ({g1} - 1) = {sv1}")
print(f"||Σ_{g2}|| = 4 * ({g2} - 1) = {sv2}")
print(f"C({total_dim}, {m}) = {coeff}")
print("\nPutting it all together:")
print(f"||Σ_{g1} x Σ_{g2}|| = {coeff} * {sv1} * {sv2} = {final_sv}")
