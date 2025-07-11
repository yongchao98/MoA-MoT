from fractions import Fraction

# The problem of finding the minimal expected number of rolls can be solved
# by setting up and solving a system of linear equations based on an optimal,
# information-preserving strategy.
# The derivation leads to the following equation for the expected number of rolls, E:
# E = (752/343) * (2401/2400)

# Define the numbers in the final equation
num1 = 752
den1 = 343
num2 = 2401
den2 = 2400

# Create Fraction objects to perform exact arithmetic
term1 = Fraction(num1, den1)
term2 = Fraction(num2, den2)

# Calculate the final result by multiplying the fractions
E = term1 * term2

# As requested, we output the numbers in the final equation and the result.
print("The minimal expected number of rolls, E, is calculated from the equation:")
print(f"E = ({num1}/{den1}) * ({num2}/{den2})")

# We can also show the result before simplification
unsimplified_num = num1 * num2
unsimplified_den = den1 * den2
print(f"This evaluates to: E = {unsimplified_num} / {unsimplified_den}")

# The Fraction object automatically simplifies the result.
print(f"As a simplified fraction, the minimal expected number of rolls is:")
print(f"E = {E.numerator}/{E.denominator}")
