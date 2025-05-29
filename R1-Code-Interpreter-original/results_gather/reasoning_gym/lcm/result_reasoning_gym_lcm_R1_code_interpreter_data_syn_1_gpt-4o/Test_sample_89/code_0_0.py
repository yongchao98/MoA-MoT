import math

# Define the numbers
a = 43
b = 37

# Calculate the LCM using the formula
lcm = abs(a * b) // math.gcd(a, b)

# Print the result
print(lcm)