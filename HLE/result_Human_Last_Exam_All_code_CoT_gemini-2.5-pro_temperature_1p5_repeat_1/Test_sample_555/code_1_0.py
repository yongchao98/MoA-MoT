from fractions import Fraction

# This script calculates the minimal expected number of rolls for the described problem.
# The derivation leads to the equation: E = (752/343) * (2401/2400),
# where E is the minimal expected number of rolls.

# The numbers from the final equation for E
num1, den1 = 752, 343
num2, den2 = 2401, 2400

# We use the Fraction class from Python's standard library to perform the calculation
# and get the simplified result automatically.
fraction1 = Fraction(num1, den1)
fraction2 = Fraction(num2, den2)
result_fraction = fraction1 * fraction2

print("The minimal expected number of rolls, E, can be found by solving a system of recurrence relations.")
print("The final equation for E is:")
# The user was asked to output each number in the final equation.
print(f"E = ({num1}/{den1}) * ({num2}/{den2})")

# Calculate and show the unsimplified fraction
unsimplified_num = num1 * num2
unsimplified_den = den1 * den2
print(f"E = {unsimplified_num}/{unsimplified_den}")

# Display the final, simplified result
final_num = result_fraction.numerator
final_den = result_fraction.denominator
print(f"The result as a simplified fraction is: {final_num}/{final_den}")