import math

# You can change this value to calculate the risk for a different 'n'
n = 10

print(f"The minimax risk for n={n} is given by the formula: 1 / (4 * (n + 1)^2)")

# Define the numbers that appear in the final equation
numerator = 1
four = 4
n_val = n
one = 1
power = 2

# Perform the calculation step-by-step
term_in_parentheses = n_val + one
term_squared = term_in_parentheses**power
denominator_val = four * term_squared
risk = numerator / denominator_val

# Output the equation with the numbers substituted, showing the calculation
print("\nSubstituting the values into the formula:")
print(f"Risk = {numerator} / ({four} * ({n_val} + {one})^{power})")
print(f"     = {numerator} / ({four} * ({term_in_parentheses})^{power})")
print(f"     = {numerator} / ({four} * ({term_squared}))")
print(f"     = {numerator} / {denominator_val}")
print(f"\nThe final minimax risk is: {risk}")
