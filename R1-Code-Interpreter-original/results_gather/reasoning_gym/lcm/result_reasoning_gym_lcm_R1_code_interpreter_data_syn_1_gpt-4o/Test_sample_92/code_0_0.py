import math

# Given numbers
a = 8
b = 54

# Calculate GCD
gcd = math.gcd(a, b)

# Calculate LCM using the formula
lcm = abs(a * b) // gcd

# Output the result
print(lcm)