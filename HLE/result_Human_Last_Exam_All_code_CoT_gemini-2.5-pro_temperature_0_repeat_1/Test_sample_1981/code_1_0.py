import math

# Based on the derivation and correction of a likely typo in the problem statement,
# the value of X0 is determined to be 1.
X0 = 1

# The expression to be calculated is 10^30 * X0^2 - 10^30 * X0 + 10.
# We substitute X0 = 1 into the expression.

# Define the terms of the equation
# Using float for calculation as the terms cancel out, avoiding large number arithmetic.
term1 = 10**30 * X0**2
term2 = 10**30 * X0
term3 = 10

# Calculate the final result
result = term1 - term2 + term3

# As requested, we output each number in the final equation and the result.
# We format the output to show the equation with the value of X0 substituted.
print(f"10**30 * ({X0})**2 - 10**30 * ({X0}) + 10 = {result}")
