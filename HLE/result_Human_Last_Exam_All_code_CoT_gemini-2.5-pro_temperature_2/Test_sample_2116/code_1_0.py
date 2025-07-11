import math

# The problem requires calculating the expected maximum earthquake magnitude.
# Through mathematical derivation, the final formula for this expectation (E) is found to be:
# E = pi / (2 * ln(2))
# This script calculates the numerical value of this expression.

# Define the constants from the formula
pi = math.pi
two = 2.0
ln2 = math.log(2)  # Natural logarithm of 2

# Calculate the final result
expected_magnitude = pi / (two * ln2)

print("The formula for the expected maximum earthquake magnitude is:")
print("E = pi / (2 * ln(2))")
print("\nWhere the values are:")
print(f"pi = {pi}")
print(f"2 = {two}")
print(f"ln(2) = {ln2}")
print("\nPlugging the numbers into the equation:")
print(f"E = {pi} / ({two} * {ln2})")
print("\nThe calculated expected maximum earthquake magnitude is:")
print(f"E = {expected_magnitude}")
