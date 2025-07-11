import math

# Based on the analytical solution, the expression simplifies to:
# (3/2) * 10^(10/3) + 37/4

# Define the constants from the simplified expression
c1 = 3/2
base = 10
exponent = 10/3
c2 = 37/4

# Calculate the final result
result = c1 * (base ** exponent) + c2

# Print the final equation with each number and the calculated result
print(f"The final calculation is:")
print(f"({c1}) * {base}**({exponent:.4f}) + ({c2}) = {result}")
