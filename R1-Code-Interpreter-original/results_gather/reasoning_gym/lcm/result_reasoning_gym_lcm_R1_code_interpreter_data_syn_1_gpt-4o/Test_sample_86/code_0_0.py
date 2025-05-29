import math

# Define the numbers
a = 39
b = 43

# Calculate the LCM using the formula
lcm = abs(a * b) // math.gcd(a, b)

# Print the result
print(lcm)