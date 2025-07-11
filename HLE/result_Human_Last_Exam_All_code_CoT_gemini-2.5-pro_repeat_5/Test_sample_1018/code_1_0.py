import math

# Step 1: Define the genera and dimensions of the surfaces.
g1 = 31
g2 = 17
m = 2  # Dimension of Sigma_31
n = 2  # Dimension of Sigma_17

# Step 2: Calculate the simplicial volume for each surface.
# The formula for the simplicial volume of a surface of genus g >= 2 is ||Sigma_g|| = 4g - 4.
sv1 = 4 * g1 - 4
sv2 = 4 * g2 - 4

# Step 3: Calculate the binomial coefficient for the product formula.
# The formula is C(m+n, m).
binom_coeff = math.comb(m + n, m)

# Step 4: Calculate the final simplicial volume of the product.
# ||Sigma_g1 x Sigma_g2|| = C(m+n, m) * ||Sigma_g1|| * ||Sigma_g2||
final_sv = binom_coeff * sv1 * sv2

# Step 5: Print the final equation with all the components.
print(f"The simplicial volume of Sigma_{g1} is ||Sigma_{g1}|| = 4 * {g1} - 4 = {sv1}")
print(f"The simplicial volume of Sigma_{g2} is ||Sigma_{g2}|| = 4 * {g2} - 4 = {sv2}")
print(f"The binomial coefficient is C({m}+{n}, {m}) = {binom_coeff}")
print(f"The simplicial volume of the product is ||Sigma_{g1} x Sigma_{g2}|| = {binom_coeff} * {sv1} * {sv2} = {final_sv}")