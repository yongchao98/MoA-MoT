import math

# Given numbers
numerator = 3375
denominator = 8400

# Calculate GCD
gcd = math.gcd(numerator, denominator)

# Simplify the fraction
simplified_numerator = numerator // gcd
simplified_denominator = denominator // gcd

# Output the simplified fraction
print(f"{simplified_numerator}/{simplified_denominator}")