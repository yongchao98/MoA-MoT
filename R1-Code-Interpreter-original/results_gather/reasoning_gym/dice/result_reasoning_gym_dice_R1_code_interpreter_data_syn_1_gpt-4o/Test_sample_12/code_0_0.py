from itertools import product
from math import gcd

# Define the number of sides for each die
sides = [20, 19, 18, 15]

# Calculate total possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate favorable outcomes
favorable_outcomes = 0
for roll in product(range(1, 21), range(1, 20), range(1, 19), range(1, 16)):
    if sum(roll) >= 48:
        favorable_outcomes += 1

# Calculate the probability
numerator = favorable_outcomes
denominator = total_outcomes

# Reduce the fraction
def reduce_fraction(n, d):
    common_divisor = gcd(n, d)
    return n // common_divisor, d // common_divisor

reduced_numerator, reduced_denominator = reduce_fraction(numerator, denominator)

# Output the reduced fraction
print(f"{reduced_numerator}/{reduced_denominator}")