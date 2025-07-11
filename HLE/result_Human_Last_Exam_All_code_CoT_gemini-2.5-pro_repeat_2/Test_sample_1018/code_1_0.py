import math

# Define the genera of the two surfaces
g1 = 31
g2 = 17

# The dimensions of the surfaces are 2
n1 = 2
n2 = 2

# Step 1: Calculate the simplicial volume of each surface
# The formula for a surface with genus g >= 2 is ||Σ_g|| = 4 * (g - 1)
sv1 = 4 * (g1 - 1)
sv2 = 4 * (g2 - 1)

print(f"The simplicial volume of the surface Σ_g with genus g >= 2 is ||Σ_g|| = 4 * (g - 1).")
print(f"For Σ_{g1}, the simplicial volume is ||Σ_{g1}|| = 4 * ({g1} - 1) = {sv1}.")
print(f"For Σ_{g2}, the simplicial volume is ||Σ_{g2}|| = 4 * ({g2} - 1) = {sv2}.")
print("-" * 20)

# Step 2: Calculate the binomial coefficient for the product formula
# The formula is C(dim(M1) + dim(M2), dim(M1))
coefficient = math.comb(n1 + n2, n1)

print(f"The product formula for simplicial volume is ||Σ_{g1} x Σ_{g2}|| = C({n1} + {n2}, {n1}) * ||Σ_{g1}|| * ||Σ_{g2}||.")
print(f"The binomial coefficient C({n1 + n2}, {n1}) is {coefficient}.")
print("-" * 20)

# Step 3: Calculate the final simplicial volume of the product space
result = coefficient * sv1 * sv2

print("The final calculation is:")
print(f"||\u03A3_{g1} x \u03A3_{g2}|| = {coefficient} * {sv1} * {sv2} = {result}")
