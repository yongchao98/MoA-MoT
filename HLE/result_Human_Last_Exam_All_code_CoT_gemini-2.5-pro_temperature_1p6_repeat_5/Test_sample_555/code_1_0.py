from fractions import Fraction

# We solve the system of linear equations derived from the optimal strategy.
# Let E1 be the minimal expected number of rolls.
# The system is:
# E1 = 1 + E7
# E7 = 1 + (9/49) * E9
# E9 = 1 + (1/21) * E3
# E3 = 1 + (1/21) * E1

# We use substitution to solve for E1.
# First, express E3 in terms of E1. Let E_M = c_M + d_M * E1.
# From E3 = 1 + (1/21) * E1:
c3 = Fraction(1)
d3 = Fraction(1, 21)

# Substitute into the equation for E9:
# E9 = 1 + (1/21) * E3 = 1 + (1/21) * (c3 + d3 * E1)
c9 = 1 + Fraction(1, 21) * c3
d9 = Fraction(1, 21) * d3

# Substitute into the equation for E7:
# E7 = 1 + (9/49) * E9 = 1 + (9/49) * (c9 + d9 * E1)
c7 = 1 + Fraction(9, 49) * c9
d7 = Fraction(9, 49) * d9

# Finally, substitute into the equation for E1:
# E1 = 1 + E7 = 1 + (c7 + d7 * E1)
# E1 * (1 - d7) = 1 + c7
# E1 = (1 + c7) / (1 - d7)

# Now, we calculate the final values for the numerator and denominator.
numerator = 1 + c7
denominator = 1 - d7

# Calculate the unsimplified numerator and denominator of the final fraction for E1
unsimplified_numerator = numerator.numerator * denominator.denominator
unsimplified_denominator = numerator.denominator * denominator.numerator

# The final result as a simplified fraction
result = numerator / denominator

print("The minimal expected number of rolls, E, is found by solving a system of linear equations.")
print("The final equation for E, derived from substitutions, is:")
print(f"E = (1 + {c7.numerator}/{c7.denominator}) / (1 - {d7.numerator}/{d7.denominator})")
print(f"E = ({numerator.numerator}/{numerator.denominator}) / ({denominator.numerator}/{denominator.denominator})")
print(f"E = {unsimplified_numerator} / {unsimplified_denominator}")
print("\nAs a simplified fraction, the result is:")
print(f"E = {result.numerator} / {result.denominator}")
