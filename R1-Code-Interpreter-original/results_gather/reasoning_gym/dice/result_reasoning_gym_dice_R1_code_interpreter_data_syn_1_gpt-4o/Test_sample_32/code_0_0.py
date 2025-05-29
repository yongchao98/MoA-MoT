from itertools import product

# Define the number of sides for each die
sides = [20, 17, 14, 5]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, sides[0] + 1), range(1, sides[1] + 1), range(1, sides[2] + 1), range(1, sides[3] + 1)):
    if sum(roll) >= 38:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from math import gcd

numerator = favorable_outcomes
denominator = total_outcomes
common_divisor = gcd(numerator, denominator)

# Reduce the fraction
reduced_numerator = numerator // common_divisor
reduced_denominator = denominator // common_divisor

# Print the result
print(f"{reduced_numerator}/{reduced_denominator}")