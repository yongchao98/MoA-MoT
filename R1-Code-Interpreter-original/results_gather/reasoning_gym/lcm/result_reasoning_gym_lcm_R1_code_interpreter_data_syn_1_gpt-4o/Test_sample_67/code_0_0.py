import math

# Define the numbers
a = 30
b = 32

# Calculate the GCD
gcd = math.gcd(a, b)

# Calculate the LCM using the formula
lcm = abs(a * b) // gcd

# Print the result
print(lcm)