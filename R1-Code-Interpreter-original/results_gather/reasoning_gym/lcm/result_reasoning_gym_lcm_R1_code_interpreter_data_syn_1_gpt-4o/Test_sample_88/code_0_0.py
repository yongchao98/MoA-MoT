import math

# Define the numbers
a = 29
b = 33

# Calculate the LCM using the formula
lcm = abs(a * b) // math.gcd(a, b)

# Print the result
print(lcm)