import math

# This code calculates the numerical coefficients for the total cross section formula.

# The formula for the total cross section (sigma) in the high-energy limit is:
# sigma = (1/2) * ( <|M|^2> / (64 * pi^2 * s) ) * (4 * pi)

# In the high-energy limit:
# <|M|^2> = 3 * g^4
# s = 4 * E^2
# g is the coupling constant, E is the center-of-mass energy of one fermion.

# Substitute these into the formula for sigma:
# sigma = (1/2) * ( (3 * g^4) / (64 * pi^2 * (4 * E^2)) ) * (4 * pi)
# sigma = (1/2) * (3 * g^4) / (256 * pi^2 * E^2) * (4 * pi)

# Let's determine the numerical coefficients.
numerator_coeff = 1 * 3 * 4
denominator_coeff = 2 * 256

# Simplify the fraction by dividing by the greatest common divisor.
common_divisor = math.gcd(numerator_coeff, denominator_coeff)
final_num = numerator_coeff // common_divisor
final_den = denominator_coeff // common_divisor

# Now we print the final formula, showing each number explicitly.
print("The final expression for the total cross section sigma in the high-energy limit is:")
print(f"sigma = ({final_num} * g**4) / ({final_den} * pi * E**2)")