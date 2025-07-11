import math

# Define the given values for n and m.
n = 54669191328000
m = 1470432000

# As derived in the thinking steps, the groupoid cardinality of the inserter
# Ins(F,G) is given by the fraction m/n.
# We need to simplify this fraction to its lowest terms.

# Calculate the greatest common divisor (GCD) of m and n.
common_divisor = math.gcd(m, n)

# Divide both the numerator and the denominator by their GCD to simplify.
numerator = m // common_divisor
denominator = n // common_divisor

# The problem asks to output the numbers in the final equation.
# So we print the simplified fraction.
print(f"{numerator} / {denominator}")
