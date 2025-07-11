import math

# The problem simplifies to calculating the value of the expression:
# (3/2) * 10^(10/3) + 37/4

# Define the numbers in the final equation
c1 = 3/2
base = 10
exponent = 10/3
c2 = 37/4

# Calculate the result
result = c1 * (base**exponent) + c2

# Output the final equation with all its numerical components and the result
print(f"The final simplified equation is:")
print(f"({c1}) * {base}^({exponent}) + ({c2})")
print(f"The result is: {result}")