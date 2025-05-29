import math

# Given numbers
a = 15
b = 12

# Calculate GCD
gcd = math.gcd(a, b)

# Calculate LCM using the formula
lcm = abs(a * b) // gcd

# Output the result
print(lcm)