import math

# The value of the integral I is pi * ln(1 + sqrt(2)).
# This script calculates the numerical value of this expression,
# showing the components as requested.

# The final equation is I = pi * ln(1 + sqrt(2)).
# The numbers in this equation are:
pi_val = math.pi
one = 1.0
two = 2.0

# Calculate the result step-by-step
sqrt_of_two = math.sqrt(two)
argument_of_ln = one + sqrt_of_two
ln_value = math.log(argument_of_ln)
final_result = pi_val * ln_value

# Output the equation and its components
print(f"The integral evaluates to the expression: pi * ln(1 + sqrt(2))")
print(f"Let's print the numbers in the final equation:")
print(f"pi = {pi_val}")
print(f"1 = {one}")
print(f"2 = {two}")
print("\nNow, we compute the final numerical value:")
print(f"I = {pi_val} * ln({one} + {sqrt_of_two})")
print(f"I = {final_result}")
