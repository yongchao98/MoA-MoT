from itertools import product

# Define the range of each die
d20 = range(1, 21)
d8 = range(1, 9)

# Total number of possible outcomes
total_outcomes = 20 * 8 * 8 * 8

# Count favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll in product(d20, d8, d8, d8):
    if sum(roll) >= 17:
        favorable_outcomes += 1

# Calculate probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability
print(probability)