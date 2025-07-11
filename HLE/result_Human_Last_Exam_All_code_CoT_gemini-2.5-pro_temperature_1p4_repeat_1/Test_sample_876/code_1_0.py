import math

# --- Given parameters and constants ---
r_w = 11.25  # Radius of the wire in mm
m = 17.32    # Distance from wire center to cable center in mm
R = 32.32    # External radius of the cable (inner radius of screen) in mm
epsilon_r = 4.2  # Relative permittivity of the insulator
epsilon_0 = 8.854e-12  # Permittivity of free space in F/m

# --- Calculation ---
# The formula for capacitance per phase is:
# C = (2 * pi * epsilon_0 * epsilon_r) / ln( (R/r_w) * (R^2 - m^2) / (R^2 + 2*m^2) )

# Calculate the argument of the natural logarithm
# Units (mm) cancel out, so we can use them directly here
log_argument = (R / r_w) * (R**2 - m**2) / (R**2 + 2 * m**2)

# Calculate the capacitance in Farads per meter (F/m)
capacitance_F_per_m = (2 * math.pi * epsilon_0 * epsilon_r) / math.log(log_argument)

# Convert from F/m to microfarads per kilometer (uF/km)
# 1 F/m = 10^9 uF/km
capacitance_uF_per_km = capacitance_F_per_m * 1e9

# --- Output the results ---
print("This script calculates the capacitance of a three-phase cable with a common screen.")
print("The formula used is C = (2 * pi * epsilon_0 * epsilon_r) / ln( (R/r_w) * (R^2 - m^2) / (R^2 + 2*m^2) )")
print("\n--- Final Equation with Substituted Values ---")

# Using f-string to embed the numerical values directly into the equation string
equation_str = (
    f"C (uF/km) = (2 * {math.pi:.5f} * {epsilon_0:.3e} F/m * {epsilon_r}) / "
    f"ln( ({R}/{r_w}) * ({R**2:.2f} - {m**2:.2f}) / ({R**2:.2f} + 2*{m**2:.2f}) ) * 1e9"
)
print(equation_str)

print("\n--- Calculation Result ---")
print(f"The capacitance of the cable is: {capacitance_uF_per_km:.3f} uF/km")
