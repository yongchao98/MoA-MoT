import math

# Step 1: Define the genera and dimensions of the surfaces.
g1 = 31
g2 = 17
m = 2
n = 2

print(f"We are computing the simplicial volume of Σ_{g1} x Σ_{g2}.")
print("-" * 20)

# Step 2: Compute the simplicial volume of the first surface, Σ_g1.
# The formula for the simplicial volume of Σ_g (for g >= 2) is ||Σ_g|| = 2 * (2g - 2).
sv1 = 2 * (2 * g1 - 2)
print(f"For the first surface, Σ_{g1}:")
print(f"||Σ_{g1}|| = 2 * (2 * {g1} - 2) = {sv1}")
print("-" * 20)

# Step 3: Compute the simplicial volume of the second surface, Σ_g2.
sv2 = 2 * (2 * g2 - 2)
print(f"For the second surface, Σ_{g2}:")
print(f"||Σ_{g2}|| = 2 * (2 * {g2} - 2) = {sv2}")
print("-" * 20)

# Step 4: Use the product formula for simplicial volumes.
# The formula is ||M x N|| = C(m+n, m) * ||M|| * ||N||.
# First, calculate the binomial coefficient C(m+n, m).
coeff = math.comb(m + n, m)
print("Using the product formula: ||M x N|| = C(m+n, m) * ||M|| * ||N||")
print(f"The dimensions are m = {m} and n = {n}.")
print(f"The coefficient C({m+n}, {m}) = {coeff}")
print("-" * 20)

# Step 5: Calculate the final simplicial volume.
total_sv = coeff * sv1 * sv2
print("The final calculation is:")
print(f"||Σ_{g1} x Σ_{g2}|| = {coeff} * {sv1} * {sv2}")
print(f"Result = {total_sv}")
