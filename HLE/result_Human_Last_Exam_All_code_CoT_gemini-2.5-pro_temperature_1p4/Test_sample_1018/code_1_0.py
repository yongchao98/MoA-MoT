import math

# Define the genera of the two surfaces
g1 = 31
g2 = 17

# The formula for the simplicial volume of a hyperbolic surface Sigma_g is ||Sigma_g|| = 4g - 4.
# We calculate this for each surface.
# Both g1 and g2 are >= 2, so the formula applies.
sv1 = 4 * g1 - 4
sv2 = 4 * g2 - 4

# For the product of two surfaces (each of dimension 2), the Gromov product formula is:
# ||Sigma_g1 x Sigma_g2|| = C(4, 2) * ||Sigma_g1|| * ||Sigma_g2||
# Calculate the binomial coefficient C(4, 2).
coeff = math.comb(4, 2)

# Calculate the final simplicial volume of the product manifold.
total_sv = coeff * sv1 * sv2

# Print the final equation with all the numbers.
print("The simplicial volume of Σ_31 x Σ_17 is calculated as C(4, 2) * ||Σ_31|| * ||Σ_17||")
print(f"||Σ_31|| = 4 * {g1} - 4 = {sv1}")
print(f"||Σ_17|| = 4 * {g2} - 4 = {sv2}")
print(f"C(4, 2) = {coeff}")
print("Final Calculation:")
print(f"{coeff} * {sv1} * {sv2} = {total_sv}")