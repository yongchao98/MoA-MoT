from itertools import product

# Define the number of sides for each die
sides = [20, 18, 14, 9]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Count the number of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, sides[0] + 1), range(1, sides[1] + 1), range(1, sides[2] + 1), range(1, sides[3] + 1)):
    if sum(roll) >= 34:
        successful_outcomes += 1

# Calculate the probability as a reduced fraction
from math import gcd

numerator = successful_outcomes
denominator = total_outcomes
common_divisor = gcd(numerator, denominator)

# Reduce the fraction
reduced_numerator = numerator // common_divisor
reduced_denominator = denominator // common_divisor

print(f"{reduced_numerator}/{reduced_denominator}")