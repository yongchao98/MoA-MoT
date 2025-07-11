# Define the genera of the two surfaces
g1 = 31
g2 = 17

# Calculate the simplicial volume for the first surface, Sigma_g1
# The formula for the simplicial volume of a surface of genus g >= 2 is 4*(g-1)
sv1 = 4 * (g1 - 1)

# Calculate the simplicial volume for the second surface, Sigma_g2
sv2 = 4 * (g2 - 1)

# The simplicial volume of the product of two surfaces is the product of their individual simplicial volumes
total_sv = sv1 * sv2

# Print the calculation steps and the final result
print(f"The simplicial volume of a surface of genus g >= 2 is given by the formula ||Σ_g|| = 4 * (g - 1).")
print(f"First, we compute the simplicial volume of Σ_{g1}:")
print(f"||Σ_{g1}|| = 4 * ({g1} - 1) = {sv1}")
print(f"Next, we compute the simplicial volume of Σ_{g2}:")
print(f"||Σ_{g2}|| = 4 * ({g2} - 1) = {sv2}")
print(f"The simplicial volume of the product Σ_{g1} x Σ_{g2} is the product of their individual volumes:")
print(f"||Σ_{g1} x Σ_{g2}|| = ||Σ_{g1}|| * ||Σ_{g2}|| = {sv1} * {sv2} = {total_sv}")