import math

# Given numbers
a = 48
b = 4

# Calculate GCD
gcd = math.gcd(a, b)

# Calculate LCM using the formula
lcm = abs(a * b) // gcd

# Output the result
print(lcm)