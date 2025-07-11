from fractions import Fraction

# This script demonstrates the calculation the parrot can perform.

# 1. Define the components of the calculation as fractions
# All numerators and denominators must be integers <= 10.
vol_const = Fraction(4, 3)
pi_approx = Fraction(10, 3)  # Approximation for pi
radius = Fraction(1, 2)
density = Fraction(9, 10)

# 2. Perform the calculation for the mass
# M = (4/3) * pi * r^3 * density
estimated_mass = vol_const * pi_approx * (radius**3) * density

# 3. Print the calculation in a way the parrot would see it.
# The final code needs to output each number in the final equation.
print("The calculation to estimate the mass is:")
print(f"Mass = ({vol_const.numerator}/{vol_const.denominator}) * ({pi_approx.numerator}/{pi_approx.denominator}) * ({radius.numerator}/{radius.denominator}) * ({radius.numerator}/{radius.denominator}) * ({radius.numerator}/{radius.denominator}) * ({density.numerator}/{density.denominator})")

# 4. Print the final result as a fraction.
print(f"\nThe estimated mass is: {estimated_mass.numerator}/{estimated_mass.denominator} kg")