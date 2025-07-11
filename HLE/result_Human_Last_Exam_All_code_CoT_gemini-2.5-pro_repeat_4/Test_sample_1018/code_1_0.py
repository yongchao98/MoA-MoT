import math

# Step 1: Define the genera and dimensions
g1 = 31
g2 = 17
m = 2
n = 2

# Step 2: Calculate the simplicial volume of each surface
# The formula is ||Σ_g|| = 4g - 4 for g > 1
sv1 = 4 * g1 - 4
sv2 = 4 * g2 - 4

print(f"The simplicial volume of a surface Σ_g is ||Σ_g|| = 4g - 4 (for g > 1).")
print(f"For Σ_{g1}, the simplicial volume is 4 * {g1} - 4 = {sv1}.")
print(f"For Σ_{g2}, the simplicial volume is 4 * {g2} - 4 = {sv2}.")
print("-" * 20)

# Step 3: Calculate the binomial coefficient for the product formula
# The formula is ||M x N|| = C(m+n, m) * ||M|| * ||N||
c = math.comb(m + n, m)

print(f"The product formula for simplicial volume is ||Σ_{g1} x Σ_{g2}|| = C({m+n}, {m}) * ||Σ_{g1}|| * ||Σ_{g2}||.")
print(f"The binomial coefficient C({m+n}, {m}) is {c}.")
print("-" * 20)

# Step 4: Calculate the final simplicial volume of the product
result = c * sv1 * sv2

# Step 5: Print the final calculation and result
print("The final calculation is:")
print(f"||Σ_{g1} x Σ_{g2}|| = {c} * {sv1} * {sv2} = {result}")
