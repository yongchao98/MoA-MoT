import math

# Given initial values
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# The formula for energy loss per centimeter (-dE/dx) is derived from Geiger's rule:
# -dE/dx = (2 * E0) / (3 * R0^(2/3) * (R0 - x)^(1/3))

# Print the formula and the substituted values
print("The formula for energy loss (-dE/dx) is:")
print("-dE/dx = (2 * E0) / (3 * R0^(2/3) * (R0 - x)^(1/3))")
print("\nSubstituting the given values:")
print(f"E0 = {E0} MeV")
print(f"R0 = {R0} cm")
print(f"x = {x} cm")
print("\nThe equation becomes:")
print(f"-dE/dx = (2 * {E0}) / (3 * {R0}^(2/3) * ({R0} - {x})^(1/3))")

# Calculate the terms in the equation
numerator = 2 * E0
r0_pow_2_3 = R0**(2/3)
remaining_range = R0 - x
remaining_range_pow_1_3 = remaining_range**(1/3)
denominator = 3 * r0_pow_2_3 * remaining_range_pow_1_3

# Perform the final calculation
energy_loss_per_cm = numerator / denominator

print("\nStep-by-step calculation:")
print(f"Numerator = 2 * {E0} = {numerator}")
print(f"Term in denominator 1: R0^(2/3) = {R0}^(2/3) = {r0_pow_2_3:.4f}")
print(f"Term in denominator 2: (R0 - x)^(1/3) = ({remaining_range})^(1/3) = {remaining_range_pow_1_3:.4f}")
print(f"Denominator = 3 * {r0_pow_2_3:.4f} * {remaining_range_pow_1_3:.4f} = {denominator:.4f}")

print(f"\nFinal Result: -dE/dx = {numerator} / {denominator:.4f}")
print(f"The energy loss per centimetre is {energy_loss_per_cm:.2f} MeV/cm.")