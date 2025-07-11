import math

# Step 1: Define the genera of the surfaces.
g1 = 31
g2 = 17

# Step 2 & 3: Define the formula for the simplicial volume of a surface
# ||Σ_g|| = 4g - 4 for g >= 2, and calculate the volumes.
vol_g1 = 4 * g1 - 4
vol_g2 = 4 * g2 - 4

# Step 4 & 5: Define the formula for the product and compute the final volume.
# ||Σ_g1 x Σ_g2|| = (1/6) * ||Σ_g1|| * ||Σ_g2||
# We use integer arithmetic as the result is an integer.
final_volume = (vol_g1 * vol_g2) // 6

# Print the calculation step-by-step
print(f"The simplicial volume of a surface Σ_g with genus g≥2 is ||Σ_g|| = 4g - 4.")
print(f"For g1 = {g1}, the simplicial volume is ||Σ_{g1}|| = 4 * {g1} - 4 = {vol_g1}.")
print(f"For g2 = {g2}, the simplicial volume is ||Σ_{g2}|| = 4 * {g2} - 4 = {vol_g2}.")
print(f"The simplicial volume of the product is given by the formula: ||Σ_{g1} x Σ_{g2}|| = (1/6) * ||Σ_{g1}|| * ||Σ_{g2}||.")
print(f"Plugging in the numbers, we get the final equation:")
print(f"||Σ_{g1} x Σ_{g2}|| = (1/6) * {vol_g1} * {vol_g2} = {final_volume}")