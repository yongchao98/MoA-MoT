import math

# Given values
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Calculate the remaining range
R_rem = R0 - x

# Calculate the energy loss per centimeter using the derived formula
# Formula: -dE/dx = (2 * E0 / (3 * R0)) * (R0 / (R0 - x))^(1/3)
energy_loss = (2 * E0 / (3 * R0)) * (R0 / R_rem)**(1/3)

# Print the explanation and the final result
print("This script calculates the energy loss per centimetre for an α-particle in standard air.")
print("The calculation is based on Geiger's rule for the range-energy relationship.")
print("\nGiven values:")
print(f"Initial energy (E₀): {E0} MeV")
print(f"Total range (R₀): {R0} cm")
print(f"Distance from source (x): {x} cm")

print("\nThe formula for energy loss (-dE/dx) at a distance x is:")
print("-dE/dx = (2 * E₀ / (3 * R₀)) * (R₀ / (R₀ - x))^(1/3)")

print("\nSubstituting the given values into the formula:")
print(f"-dE/dx = (2 * {E0} / (3 * {R0})) * ({R0} / ({R0} - {x}))**(1/3)")

print(f"\nResult:")
print(f"The energy loss per centimetre at a distance of {x} cm is: {energy_loss:.4f} MeV/cm")

# The final numerical answer in the required format will be the rounded value of energy_loss
# to two decimal places.
# energy_loss is approx 0.8502
# Rounded to two decimal places: 0.85
# So the final answer is 0.85
final_answer = round(energy_loss, 2)