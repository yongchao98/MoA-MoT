import math

# Based on the analysis of the problem, a likely typo exists in the definition of beta.
# Assuming the typo is corrected, the value of X0 is determined to be 1.
# We now calculate the final expression with X0 = 1.

# The expression to calculate is: 10^30 * X0^2 - 10^30 * X0 + 10
X0 = 1

# Define the coefficients
a = 10**30
b = 10**30
c = 10

# Calculate the result
result = a * (X0**2) - b * X0 + c

# Print the final equation with all its components as requested.
# Using scientific notation for large numbers to keep the output readable.
print(f"{a:.0e} * {X0}^2 - {b:.0e} * {X0} + {c} = {int(result)}")
