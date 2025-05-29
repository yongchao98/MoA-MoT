import math

# Given numbers
a = 27
b = 35

# Calculate GCD
gcd = math.gcd(a, b)

# Calculate LCM using the relationship between GCD and LCM
lcm = abs(a * b) // gcd

# Output the LCM
print(lcm)