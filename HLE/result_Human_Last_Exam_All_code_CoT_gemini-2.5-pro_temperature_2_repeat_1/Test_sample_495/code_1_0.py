import math

# Problem: Calculate the energy loss per centimetre for α-particles.
#
# Given parameters:
# Initial energy (E0) = 8.5 MeV
# Range in standard air (R) = 8.3 cm
# Distance from the source (x) = 4.0 cm

# Define the variables
E0 = 8.5
R = 8.3
x = 4.0

# Based on Geiger's rule (R = k*E^(3/2)), the energy loss per unit distance (-dE/dx)
# can be derived as: -dE/dx = (2 * E0) / (3 * R) * (R / (R - x))**(1/3)

# Calculate the energy loss per centimeter
energy_loss_per_cm = (2 * E0) / (3 * R) * (R / (R - x))**(1/3)

# Print the calculation steps and the final result
print("To calculate the energy loss per centimetre, we use the formula:")
print("(-dE/dx) = (2 * E0) / (3 * R) * (R / (R - x))^(1/3)\n")
print("Substituting the given values:")
# The f-string below displays the full equation with the numbers plugged in
print(f"(-dE/dx) = (2 * {E0}) / (3 * {R}) * ({R} / ({R} - {x}))^(1/3)")
print(f"\nCalculated energy loss per centimetre: {energy_loss_per_cm:.3f} MeV/cm")
