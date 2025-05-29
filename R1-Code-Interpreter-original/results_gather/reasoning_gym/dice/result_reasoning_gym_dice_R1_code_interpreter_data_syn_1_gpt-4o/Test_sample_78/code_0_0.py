from itertools import product

# Define the number of sides for each die
sides = [20, 15, 14, 9]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for outcome in product(range(1, sides[0] + 1), range(1, sides[1] + 1), range(1, sides[2] + 1), range(1, sides[3] + 1)):
    if sum(outcome) >= 37:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from math import gcd

def reduce_fraction(numerator, denominator):
    common_divisor = gcd(numerator, denominator)
    return numerator // common_divisor, denominator // common_divisor

numerator, denominator = reduce_fraction(favorable_outcomes, total_outcomes)

# Print the result as a reduced fraction
print(f"{numerator}/{denominator}")