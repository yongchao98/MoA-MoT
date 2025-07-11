import math

# The problem is to evaluate the sum S = sum_{n=2 to inf} (-1)^n * S_n / n.
# Based on the derivation explained above, the closed-form solution is:
# S = (ln(2))^2 / 2 + 1 - (pi)^2 / 12

# Define the constants and numbers that appear in the final formula.
ln2 = math.log(2)
pi = math.pi
num_1 = 1.0
num_2 = 2.0
num_12 = 12.0

# Calculate the final value of the expression.
result = (ln2**2 / num_2) + num_1 - (pi**2 / num_12)

# Print the formula and its components.
print(f"The sum evaluates to the expression: (ln(2))^2 / 2 + 1 - (pi)^2 / 12")
print(f"\nUsing the following numerical values for the components:")
print(f"ln(2) = {ln2}")
print(f"pi = {pi}")
print(f"1 = {num_1}")
print(f"2 = {num_2}")
print(f"12 = {num_12}")

# Display the final equation with the numbers plugged in, showing the result.
print("\nThe final equation with the numerical values is:")
# This line prints each number in the final equation as requested.
print(f"(({ln2})^2 / {num_2}) + {num_1} - (({pi})^2 / {num_12})")
print(f"\nWhich evaluates to the final result:")
print(f"{result}")
