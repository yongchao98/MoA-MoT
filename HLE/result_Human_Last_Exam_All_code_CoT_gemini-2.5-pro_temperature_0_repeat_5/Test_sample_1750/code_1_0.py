import math

# The problem simplifies to calculating the fraction 35 * C(50, 25) / 2^50.
# We compute the numerator and denominator and then simplify the fraction.

n = 50
k = 25

# Calculate the numerator and denominator
numerator_val = 35 * math.comb(n, k)
denominator_val = 2**n

# To simplify the fraction, we find the greatest common divisor (GCD).
common_divisor = math.gcd(numerator_val, denominator_val)

# Simplify the fraction
simplified_numerator = numerator_val // common_divisor
simplified_denominator = denominator_val // common_divisor

# Output the numbers in the final equation
print(f"The integral evaluates to (35 * C(50, 25)) / 2^50 * pi")
print(f"We are asked for the fractional part, which is (35 * {math.comb(n, k)}) / {2**n}")
print(f"Numerator = {numerator_val}")
print(f"Denominator = {denominator_val}")
print(f"After simplification, the fraction is:")
print(f"Numerator: {simplified_numerator}")
print(f"Denominator: {simplified_denominator}")
print(f"Final Answer as a fraction: {simplified_numerator}/{simplified_denominator}")
