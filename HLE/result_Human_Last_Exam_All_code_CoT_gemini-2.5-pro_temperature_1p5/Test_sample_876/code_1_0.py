import math

# Plan:
# 1. Define the given parameters for the three-phase cable.
# 2. Use the standard formula for capacitance per phase.
# 3. Calculate the capacitance in Farads per meter (F/m).
# 4. Convert the result to microfarads per kilometer (uF/km).
# 5. Print the breakdown of the calculation and the final answer.

# --- Given parameters ---
# r_w: radius of the wire in mm
r_w = 11.25
# m: distance from the center of a wire to the center of the cable in mm
m = 17.32
# R: inner radius of the common screen in mm
R = 32.32
# epsilon_r: relative permittivity of the insulator
epsilon_r = 4.2

# --- Physical Constant ---
# epsilon_0: permittivity of free space in F/m
epsilon_0 = 8.854e-12

# --- Step-by-step Calculation ---
# 1. Calculate the permittivity of the insulator
epsilon = epsilon_r * epsilon_0

# 2. Calculate the term inside the natural logarithm, ln(x).
#    The units (mm) cancel out, so no conversion is needed for this part.
x_numerator = R**2 - m**2
x_denominator = r_w * R
x = x_numerator / x_denominator

# 3. Calculate the full numerator of the capacitance formula
C_numerator = 2 * math.pi * epsilon

# 4. Calculate the capacitance in Farads per meter (F/m)
C_F_per_m = C_numerator / math.log(x)

# 5. Convert capacitance to microfarads per kilometer (uF/km)
#    Conversion factor: 1 F/m = 1e6 uF / (1e-3 km) = 1e9 uF/km
C_muF_per_km = C_F_per_m * 1e9

# --- Output the final equation and result ---
print("The capacitance (C) calculation is as follows:\n")
print("Formula: C = (2 * pi * epsilon_r * epsilon_0) / ln((R^2 - m^2) / (r_w * R))\n")
print("Substituting the given values into the formula:")
print(f"C = (2 * {math.pi:.5f} * {epsilon_r} * {epsilon_0:.3e}) / ln(({R}^2 - {m}^2) / ({r_w} * {R}))")
print(f"C = ({C_numerator:.5e}) / ln(({x_numerator:.3f}) / ({x_denominator:.3f}))")
print(f"C = {C_numerator:.5e} / ln({x:.5f})")
print(f"C = {C_numerator:.5e} / {math.log(x):.5f}")
print(f"C = {C_F_per_m:.5e} F/m")
print("\nConverting the result to microfarads per kilometer (uF/km):")
print(f"C = {C_F_per_m:.5e} F/m * 10^9 = {C_muF_per_km:.4f} uF/km")

print(f"\nFinal Result: The capacitance of the cable is {C_muF_per_km:.4f} \u03BCF/km.")
<<<0.3260>>>