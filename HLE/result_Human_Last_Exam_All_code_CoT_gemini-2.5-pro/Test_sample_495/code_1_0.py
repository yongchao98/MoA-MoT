import math

# --- Given initial conditions ---
# Initial energy in MeV
E0 = 8.5
# Total range in cm
R0 = 8.3
# Distance from the source in cm
x = 4.0

# --- Calculation ---
# The formula for energy loss per centimeter |dE/dx| at a distance x is derived from Geiger's rule (R = k*E^(3/2)).
# The formula is: |dE/dx| = (2/3) * E0 / (R0^(2/3) * (R0 - x)^(1/3))

# Calculate the denominator parts
R0_pow_2_3 = R0**(2/3)
rem_range_pow_1_3 = (R0 - x)**(1/3)
denominator = R0_pow_2_3 * rem_range_pow_1_3

# Calculate the numerator
numerator = (2/3) * E0

# Calculate the energy loss per cm
energy_loss_per_cm = numerator / denominator

# --- Output the result ---
print("This script calculates the energy loss of an alpha-particle in standard air.")
print("Based on Geiger's rule, the energy loss per centimeter |dE/dx| can be found using the formula:")
print("|dE/dx| = (2/3) * E0 / (R0^(2/3) * (R0 - x)^(1/3))\n")

print("Substituting the given values:")
print(f"Initial Energy (E0) = {E0} MeV")
print(f"Total Range (R0) = {R0} cm")
print(f"Distance (x) = {x} cm\n")

print("The equation with the numbers plugged in is:")
print(f"|dE/dx| = (2/3) * {E0} / ({R0}^(2/3) * ({R0} - {x})^(1/3))")
print(f"|dE/dx| = {numerator:.4f} / ({R0_pow_2_3:.4f} * {rem_range_pow_1_3:.4f})")
print(f"|dE/dx| = {numerator:.4f} / {denominator:.4f}\n")

print(f"The calculated energy loss per centimetre at a distance of {x} cm is: {energy_loss_per_cm:.3f} MeV/cm")

# Final answer in the specified format
# print(f"<<<{energy_loss_per_cm:.3f}>>>")