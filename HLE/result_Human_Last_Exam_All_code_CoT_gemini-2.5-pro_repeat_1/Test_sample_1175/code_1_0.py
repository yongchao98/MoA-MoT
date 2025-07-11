import math

# The expression simplifies to (3/2) * 10**(10/3) + 37/4.
# We calculate this value.

# Define the components of the simplified expression
coeff1 = 3 / 2
base = 10
exponent = 10 / 3
term2 = 37 / 4

# Calculate the final result
result = coeff1 * (base ** exponent) + term2

# Print the final equation with all its components
print(f"The calculation is: {coeff1} * {base}**({exponent}) + {term2}")
print(f"Result: {result}")