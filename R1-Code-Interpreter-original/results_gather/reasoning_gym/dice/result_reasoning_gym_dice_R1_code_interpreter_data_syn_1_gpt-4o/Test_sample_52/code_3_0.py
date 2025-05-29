def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Given fraction
numerator = 869
denominator = 1020

# Calculate GCD
common_divisor = gcd(numerator, denominator)

# Reduce the fraction
reduced_numerator = numerator // common_divisor
reduced_denominator = denominator // common_divisor

print(f"{reduced_numerator}/{reduced_denominator}")