import math

# The definite integral I evaluates to the symbolic expression: pi * ln(1 + sqrt(2))
# This code calculates the numerical value of this expression.

# The final equation is I = pi * ln(1 + sqrt(2))
# We will print the value of each number/component in this equation as requested.
final_equation_str = "I = pi * ln(1 + sqrt(2))"

# The numbers in the final equation are pi, 1, and 2.
pi_val = math.pi
one = 1
two = 2

# Calculate intermediate values
sqrt_two = math.sqrt(two)
inner_log_val = one + sqrt_two
log_val = math.log(inner_log_val)

# Calculate the final result
result = pi_val * log_val

# Print the results in a structured way
print(f"The integral evaluates to the expression: {final_equation_str}")
print("\nHere are the values of the components in this expression:")
print(f"The value of pi is: {pi_val}")
print(f"The number one is: {one}")
print(f"The number two is: {two}")
print(f"The value of sqrt(2) is: {sqrt_two}")
print(f"The argument of the natural logarithm, 1 + sqrt(2), is: {inner_log_val}")
print(f"The value of ln(1 + sqrt(2)) is: {log_val}")
print("\nFinally, the numerical value of the integral is:")
print(f"I = {result}")
