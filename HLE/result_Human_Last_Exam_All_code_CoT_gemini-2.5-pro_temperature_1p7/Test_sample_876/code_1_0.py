import math

# Define given values and constants in SI units (meters)
r_w = 11.25 / 1000  # radius of the wires in meters
m = 17.32 / 1000    # distance from wire center to cable center in meters
R = 32.32 / 1000    # external radius of the cable in meters
eps_r = 4.2         # relative permittivity of the insulator (dimensionless)
eps_0 = 8.854e-12   # permittivity of free space in F/m

# Calculate the absolute permittivity of the insulator
eps = eps_r * eps_0

# Calculate the term inside the natural logarithm
# This term must be dimensionless
log_argument = R / (r_w * (1 - (m / R)**2))

# Calculate the capacitance in Farads per meter (F/m)
capacitance_F_per_m = (2 * math.pi * eps) / math.log(log_argument)

# Convert capacitance from F/m to µF/km
# 1 F/m = 10^6 µF / (1e-3 km) = 10^9 µF/km
conversion_factor = 1e9
capacitance_uF_per_km = capacitance_F_per_m * conversion_factor

# Print the final equation with numerical values
print("The capacitance C is calculated using the formula:")
print("C (in µF/km) = (2 * pi * ε_r * ε_0) / ln[R / (r_w * (1 - (m/R)^2))] * 10^9")
print(f"C = (2 * {math.pi:.4f} * {eps_r} * {eps_0}) / ln[{R:.5f} / ({r_w:.5f} * (1 - ({m:.5f}/{R:.5f})^2))] * 10^9")

# Print the final result
print("\nFinal Result:")
print(f"The capacitance of the three-phase cable is: {capacitance_uF_per_km:.4f} µF/km")