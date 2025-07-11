import math

# Initial given values
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Calculate the remaining range
R_remaining = R0 - x

# The formula for energy loss is: (2 * E0) / (3 * R0^(2/3) * (R0 - x)^(1/3))
# Calculate the numerator and denominator separately for clarity
numerator = 2 * E0
denominator_part1 = 3
denominator_part2 = math.pow(R0, 2/3)
denominator_part3 = math.pow(R_remaining, 1/3)
denominator = denominator_part1 * denominator_part2 * denominator_part3

# Calculate the energy loss
energy_loss = numerator / denominator

# Print the final equation with all the numbers
print("Based on Geiger's rule, the energy loss |dE/dx| can be calculated with the formula:")
print("|dE/dx| = (2 * E₀) / (3 * R₀^(2/3) * (R₀ - x)^(1/3))\n")
print("Plugging in the values:")
print(f"|dE/dx| = (2 * {E0}) / (3 * ({R0})^(2/3) * ({R0} - {x})^(1/3))")
print(f"|dE/dx| = {numerator} / (3 * ({R0})^(2/3) * ({R_remaining})^(1/3))")
print(f"|dE/dx| = {numerator} / ({denominator_part1} * {denominator_part2:.4f} * {denominator_part3:.4f})")
print(f"|dE/dx| = {numerator} / {denominator:.4f}\n")

# Print the final result
print(f"The calculated energy loss per centimetre at a distance of {x} cm is: {energy_loss:.4f} MeV/cm")