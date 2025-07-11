import math

# --- Given values ---
E0 = 8.5  # Initial energy in MeV
R = 8.3   # Total range in cm
x = 4.0   # Distance from the source in cm

# --- Explanation ---
print("This script calculates the energy loss per centimeter for an alpha-particle in standard air.")
print("It uses an empirical formula derived from Geiger's rule (R = a * E^(3/2)).\n")
print("The derived formula for energy loss |dE/dx| is:")
print("|dE/dx| = (2 * E₀) / (3 * R^(2/3) * (R - x)^(1/3))\n")

print("--- Calculation ---")
print("Given values:")
print(f"Initial Energy (E₀) = {E0} MeV")
print(f"Total Range (R) = {R} cm")
print(f"Distance from source (x) = {x} cm\n")

# --- Applying the formula ---
# This section shows the formula with the actual numbers plugged in.
print("Final equation with the given values:")
print(f"|dE/dx| = (2 * {E0}) / (3 * {R}^(2/3) * ({R} - {x})^(1/3))")

# Calculate the numerical result
numerator = 2 * E0
denominator = 3 * (R**(2/3)) * ((R - x)**(1/3))
energy_loss_per_cm = numerator / denominator

print(f"|dE/dx| = {numerator:.3f} / {denominator:.3f}\n")

# --- Final Answer ---
print(f"The energy loss per centimetre at a distance of {x} cm is: {energy_loss_per_cm:.3f} MeV/cm")

<<<0.850>>>