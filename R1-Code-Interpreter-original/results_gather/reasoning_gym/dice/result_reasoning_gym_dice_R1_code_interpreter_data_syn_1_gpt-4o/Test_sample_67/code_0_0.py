from itertools import product

# Define the range of each die
d20 = range(1, 21)
d14_1 = range(1, 15)
d14_2 = range(1, 15)
d8 = range(1, 9)

# Total number of possible outcomes
total_outcomes = 20 * 14 * 14 * 8

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible combinations of dice rolls
for roll in product(d20, d14_1, d14_2, d8):
    if sum(roll) >= 26:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes).limit_denominator()

print(probability)