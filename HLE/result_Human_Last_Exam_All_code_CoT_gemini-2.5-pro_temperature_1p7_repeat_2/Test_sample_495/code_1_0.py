import math

# Given values
E0 = 8.5  # Initial energy in MeV
R = 8.3   # Total range in cm
x = 4.0   # Distance from the source in cm

# Calculate the remaining range
R_minus_x = R - x

# Calculate the energy loss per centimeter, |dE/dx|
# Formula: |dE/dx| = (2/3) * E0 * (R^(2/3))⁻¹ * ((R-x)^(1/3))⁻¹
# This can be written as: |dE/dx| = (2 * E0) / (3 * R**(2/3) * (R-x)**(1/3))
numerator = 2 * E0
denominator = 3 * (R**(2/3)) * (R_minus_x**(1/3))
energy_loss_per_cm = numerator / denominator

# Print the final equation with numbers and the result
print("The formula for energy loss per centimeter |dE/dx| is derived from Geiger's law:")
print("|dE/dx| = (2 * E₀) / (3 * R^(2/3) * (R - x)^(1/3))")
print("\nSubstituting the given values:")
print(f"|dE/dx| = (2 * {E0}) / (3 * {R}^(2/3) * ({R} - {x})^(1/3))")
print(f"|dE/dx| = {2 * E0} / (3 * {R**(2/3):.4f} * {R_minus_x**(1/3):.4f})")
print(f"|dE/dx| = {numerator} / {denominator:.4f}")
print(f"\nThe calculated energy loss per centimetre at {x} cm is: {energy_loss_per_cm:.3f} MeV/cm")
