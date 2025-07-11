import math

# Plan:
# 1. Define given parameters and constants, converting units to meters.
# 2. Use the formula for capacitance of a three-phase cable with a common screen:
#    C = (2 * pi * epsilon) / ln((R^3 - m^3) / (3 * r_w * m * R))
# 3. Calculate the capacitance in F/m.
# 4. Convert the result to µF/km.
# 5. Print the full calculation process for clarity.

# Given parameters
r_w_mm = 11.25  # radius of the wires in mm
m_mm = 17.32    # distance from the center of a wire to the center of the cable in mm
R_mm = 32.32    # external radius of the cable in mm
epsilon_r = 4.2 # relative permittivity of the insulator

# Constants
epsilon_0 = 8.854187817e-12 # Permittivity of free space in F/m

# Convert units from mm to m
r_w = r_w_mm / 1000
m = m_mm / 1000
R = R_mm / 1000

# Perform calculations
epsilon = epsilon_0 * epsilon_r
R_cubed = R**3
m_cubed = m**3
log_numerator = R_cubed - m_cubed
log_denominator = 3 * r_w * m * R
log_argument = log_numerator / log_denominator
numerator_C = 2 * math.pi * epsilon
denominator_C = math.log(log_argument)
C_per_meter = numerator_C / denominator_C
C_per_km_uF = C_per_meter * 1e9

# Print the step-by-step calculation
print("Calculation of the capacitance for a three-phase cable with a common screen.")
print("-" * 70)
print("The formula for capacitance per phase is C [F/m] = (2 * pi * epsilon_0 * epsilon_r) / ln((R^3 - m^3) / (3 * r_w * m * R))")
print("\nSubstituting the given values:")
print(f"  - r_w = {r_w_mm} mm = {r_w} m")
print(f"  - m = {m_mm} mm = {m} m")
print(f"  - R = {R_mm} mm = {R} m")
print(f"  - epsilon_r = {epsilon_r}")
print(f"  - epsilon_0 = {epsilon_0:.4e} F/m")
print(f"  - pi = {math.pi:.5f}")
print("-" * 70)

print("The full equation with numbers is:")
print(f"C [F/m] = (2 * {math.pi:.5f} * {epsilon_0:.4e} * {epsilon_r}) / ln((({R})^3 - ({m})^3) / (3 * {r_w} * {m} * {R}))")
print("\nStep-by-step evaluation:")
print(f"1. Numerator (2 * pi * epsilon):")
print(f"   2 * {math.pi:.5f} * {epsilon_0:.4e} * {epsilon_r} = {numerator_C:.4e}")
print(f"2. Denominator (ln(...)):")
print(f"   Argument = (({R:.5f})^3 - ({m:.5f})^3) / (3 * {r_w:.5f} * {m:.5f} * {R:.5f})")
print(f"            = ({R_cubed:.4e} - {m_cubed:.4e}) / ({log_denominator:.4e})")
print(f"            = {log_numerator:.4e} / {log_denominator:.4e} = {log_argument:.5f}")
print(f"   ln({log_argument:.5f}) = {denominator_C:.5f}")
print(f"\n3. Capacitance in F/m:")
print(f"   C [F/m] = {numerator_C:.4e} / {denominator_C:.5f} = {C_per_meter:.4e} F/m")
print("\n4. Final conversion to µF/km:")
print(f"   C [µF/km] = {C_per_meter:.4e} F/m * 10^9 = {C_per_km_uF:.3f} µF/km")
print("-" * 70)
print(f"The capacitance of the cable is {C_per_km_uF:.3f} µF/km.")

final_result = C_per_km_uF
print(f'<<<{final_result:.3f}>>>')