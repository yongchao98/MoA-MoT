import math

# Given numbers
a = 24
b = 22

# Calculate GCD
gcd = math.gcd(a, b)

# Calculate LCM using the formula
lcm = abs(a * b) // gcd

# Output the result
print(lcm)