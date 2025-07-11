import math

# Step 1: Define the genera of the surfaces.
g1 = 31
g2 = 17

# The dimensions of the surfaces are both 2.
n = 2
m = 2

# Step 2: Calculate the simplicial volume of each surface.
# The formula for a surface of genus g >= 2 is ||Sigma_g|| = 4g - 4.
vol1 = 4 * g1 - 4
vol2 = 4 * g2 - 4

# Step 3: Calculate the binomial coefficient for the product formula.
# The coefficient is C(n+m, n).
coefficient = math.comb(n + m, n)

# Step 4: Calculate the total simplicial volume of the product space.
total_volume = coefficient * vol1 * vol2

# Step 5: Print the final calculation, showing each number in the equation.
print("The simplicial volume of the product of two surfaces is calculated using the formula:")
print("||Sigma_g1 x Sigma_g2|| = C(n+m, n) * ||Sigma_g1|| * ||Sigma_g2||")
print(f"where ||Sigma_g|| = 4g - 4, n = dim(Sigma_g1), and m = dim(Sigma_g2).\n")
print("For this problem:")
print(f"||Sigma_{g1}|| = 4 * {g1} - 4 = {vol1}")
print(f"||Sigma_{g2}|| = 4 * {g2} - 4 = {vol2}")
print(f"C({n+m}, {n}) = {coefficient}\n")
print("The final equation is:")
print(f"||\u03A3_{g1} \u00D7 \u03A3_{g2}|| = {coefficient} \u00D7 {vol1} \u00D7 {vol2} = {total_volume}")