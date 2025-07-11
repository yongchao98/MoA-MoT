import math

# Step 1: Define the genera and dimensions
g1 = 31
g2 = 17
m = 2
n = 2

# Step 2: Calculate the simplicial volume for each surface
# The formula for the simplicial volume of a surface of genus g >= 1 is ||Σ_g|| = 4 * (g - 1).
sv1 = 4 * (g1 - 1)
sv2 = 4 * (g2 - 1)

print(f"The simplicial volume of a surface Σ_g of genus g >= 1 is given by the formula ||Σ_g|| = 4 * (g - 1).")
print(f"For the surface of genus {g1}:")
print(f"||Σ_{g1}|| = 4 * ({g1} - 1) = {sv1}")
print(f"For the surface of genus {g2}:")
print(f"||Σ_{g2}|| = 4 * ({g2} - 1) = {sv2}")
print("-" * 30)

# Step 3: Calculate the binomial coefficient for the product formula
# The product formula is ||M x N|| = C(m+n, m) * ||M|| * ||N||
# Here m = n = 2 (dimension of a surface).
binom_coeff = math.comb(m + n, m)

print(f"The product formula for simplicial volume is ||Σ_{g1} x Σ_{g2}|| = C({m+n}, {m}) * ||Σ_{g1}|| * ||Σ_{g2}||.")
print(f"The binomial coefficient C({m+n}, {m}) is {binom_coeff}.")
print("-" * 30)

# Step 4: Calculate the final simplicial volume
total_sv = binom_coeff * sv1 * sv2

# Step 5: Print the final equation with all values
print("The final calculation is:")
print(f"||Σ_{g1} x Σ_{g2}|| = {binom_coeff} * {sv1} * {sv2} = {total_sv}")
