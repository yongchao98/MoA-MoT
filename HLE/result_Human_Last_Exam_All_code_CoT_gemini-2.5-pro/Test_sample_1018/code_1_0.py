import math

# Step 1: Define the genera of the surfaces.
g1 = 31
g2 = 17

# Step 2: Define the dimensions of the surfaces.
n1 = 2
n2 = 2

# Step 3: Calculate the simplicial volume for each surface.
# The formula for the simplicial volume of a surface Sigma_g is 4g - 4.
sv1 = 4 * g1 - 4
sv2 = 4 * g2 - 4

print(f"The simplicial volume of Sigma_{g1} is ||Sigma_{g1}|| = 4 * {g1} - 4 = {sv1}")
print(f"The simplicial volume of Sigma_{g2} is ||Sigma_{g2}|| = 4 * {g2} - 4 = {sv2}")

# Step 4: Calculate the binomial coefficient for the product formula.
# The formula is C(n1+n2, n1).
C = math.comb(n1 + n2, n1)
print(f"The coefficient for the product formula is C({n1+n2}, {n1}) = {C}")

# Step 5: Calculate the simplicial volume of the product.
# ||Sigma_g1 x Sigma_g2|| = C(n1+n2, n1) * ||Sigma_g1|| * ||Sigma_g2||
result = C * sv1 * sv2

# Step 6: Print the final calculation and result.
print("\nThe final computation is:")
print(f"||Sigma_{g1} x Sigma_{g2}|| = {C} * {sv1} * {sv2} = {result}")
