from itertools import product

# Define the range of each die
d20 = range(1, 21)
d13_1 = range(1, 14)
d13_2 = range(1, 14)
d3 = range(1, 4)

# Total number of possible outcomes
total_outcomes = 20 * 13 * 13 * 3

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible combinations of rolls
for roll in product(d20, d13_1, d13_2, d3):
    if sum(roll) >= 32:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)