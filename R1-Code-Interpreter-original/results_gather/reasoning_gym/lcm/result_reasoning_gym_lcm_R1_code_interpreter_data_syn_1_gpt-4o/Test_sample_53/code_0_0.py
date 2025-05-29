import math

# Given numbers
a = 60
b = 100

# Calculate GCD
gcd = math.gcd(a, b)

# Calculate LCM using the relationship between GCD and LCM
lcm = abs(a * b) // gcd

# Output the LCM
print(lcm)