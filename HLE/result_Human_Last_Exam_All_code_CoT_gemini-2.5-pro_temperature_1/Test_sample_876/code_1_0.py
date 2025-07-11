import math

# Define the given parameters
r_w_mm = 11.25  # radius of the wire in mm
m_mm = 17.32    # distance from wire center to cable center in mm
R_mm = 32.32    # external radius of the cable screen in mm
epsilon_r = 4.2 # relative permittivity of the insulator

# Constants
epsilon_0 = 8.854e-12 # Permittivity of free space in F/m

# Convert dimensions from mm to m for calculation
r_w = r_w_mm / 1000
m = m_mm / 1000
R = R_mm / 1000

# Calculate the numerator and denominator for the natural logarithm argument
log_numerator = R**3 - m**3
log_denominator = 3 * m * r_w * R

# Calculate the argument of the natural logarithm
log_arg = log_numerator / log_denominator

# Calculate the natural logarithm term
ln_term = math.log(log_arg)

# Calculate the numerator of the capacitance formula
capacitance_numerator = 2 * math.pi * epsilon_0 * epsilon_r

# Calculate capacitance in Farads per meter (F/m)
C_F_per_m = capacitance_numerator / ln_term

# Convert capacitance from F/m to microfarads per kilometer (uF/km)
# 1 F/m = 1e6 uF / (1e-3 km) = 1e9 uF/km
C_uF_per_km = C_F_per_m * 1e9

# Print the final equation with all numbers substituted
print("The capacitance C is calculated using the formula:")
print(f"C (μF/km) = (2 * π * {epsilon_0} * {epsilon_r}) / ln((({R})**3 - ({m})**3) / (3 * {m} * {r_w} * {R})) * 1e9")
print("\nCalculation steps:")
print(f"ln argument = (({R**3:.5e}) - ({m**3:.5e})) / ({log_denominator:.5e}) = {log_arg:.4f}")
print(f"ln({log_arg:.4f}) = {ln_term:.4f}")
print(f"Capacitance in F/m = ({capacitance_numerator:.5e}) / {ln_term:.4f} = {C_F_per_m:.5e} F/m")
print("\nFinal Result:")
print(f"The capacitance of the cable is: {C_uF_per_km:.4f} μF/km")

print(f"<<<{C_uF_per_km:.4f}>>>")