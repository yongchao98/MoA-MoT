import math

# Define the numbers
a = 61
b = 61

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the formula
lcm = abs(a * b) // gcd

# Print the result
print(lcm)