import math

# --- Problem Variables ---
# E0: Initial energy of the alpha-particles in MeV
# R0: Total range of the alpha-particles in cm
# x: Distance from the source in cm
E0 = 8.5
R0 = 8.3
x = 4.0

# --- Calculation ---
# The formula for the energy loss per centimeter (-dE/dx) is derived from Geiger's Law.
# It is given by: -dE/dx = (2/3) * (E0 / R0) * ((R0 - x) / R0)^(-1/3)

# Calculate each part of the equation
constant_factor = 2/3
energy_range_ratio = E0 / R0
remaining_range = R0 - x
remaining_range_ratio = remaining_range / R0
bragg_curve_factor = remaining_range_ratio**(-1/3.0)

# Calculate the final result
energy_loss_per_cm = constant_factor * energy_range_ratio * bragg_curve_factor

# --- Output ---
print("Calculating the energy loss per centimetre (-dE/dx) for an α-particle.")
print("The formula used is: -dE/dx = (2/3) * (E0 / R0) * ((R0 - x) / R0)**(-1/3)\n")
print("Given values:")
print(f"  Initial Energy (E0) = {E0} MeV")
print(f"  Total Range (R0)    = {R0} cm")
print(f"  Distance (x)        = {x} cm\n")

print("The final equation with the numbers substituted is:")
# The f-string formats the numbers nicely within the equation string.
print(f"-dE/dx = (2/3) * ({E0} / {R0}) * (({R0} - {x}) / {R0})**(-1/3)")
print(f"-dE/dx = ({constant_factor:.4f}) * ({energy_range_ratio:.4f}) * ({bragg_curve_factor:.4f})")
print(f"-dE/dx = {energy_loss_per_cm:.4f}\n")

print(f"The calculated energy loss per centimetre at a distance of {x} cm is {energy_loss_per_cm:.3f} MeV/cm.")
