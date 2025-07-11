import math
from fractions import Fraction

# This script derives the calculation for the mass of the rock
# using 5-bit fractional arithmetic, aiming for the smallest absolute error.

# --- Problem setup ---
# Mass = Density * Volume
# Volume = 4/3 * pi * r^3
# Mass = Density * 4/3 * pi * r^3

# --- Initial Values as 5-bit Fractions ---
# Density (rho) = 0.9 kg/cm^3 -> 9/10
# Radius (r) = 0.5 cm -> 1/2
rho = Fraction(9, 10)
r = Fraction(1, 2)
four_thirds = Fraction(4, 3)

# After analysis, the best 5-bit approximation for pi for this problem is 28/9
pi_approx = Fraction(28, 9)

# --- Step-by-step calculation demonstrating the process ---
# 1. Cube the radius: r^3 = (1/2)^3 = 1/8
r_cubed = Fraction(1, 8)

# 2. Combine constant terms: rho * 4/3 * r^3
# This is done by simplifying at each step to avoid numbers > 31.
# (9/10 * 4/3) -> simplify to (6/5)
# (6/5 * 1/8) -> simplify to (3/20)
constants_part = Fraction(3, 20)

# 3. Multiply by the pi approximation: (3/20) * (28/9)
# Simplify to avoid overflow: (3/9) * (28/20) = (1/3) * (7/5) = 7/15
final_mass = Fraction(7, 15)

# --- Final Equation ---
# The full calculation, using the chosen fractions, is:
# Mass = (9/10) * (4/3) * (28/9) * (1/8)
# which simplifies down to 7/15.

print("The final calculation is represented by the equation:")
# We print each number in the equation as requested.
print(f"{rho.numerator} / {rho.denominator} * {four_thirds.numerator} / {four_thirds.denominator} * {pi_approx.numerator} / {pi_approx.denominator} * {r_cubed.numerator} / {r_cubed.denominator} = {final_mass.numerator} / {final_mass.denominator}")
