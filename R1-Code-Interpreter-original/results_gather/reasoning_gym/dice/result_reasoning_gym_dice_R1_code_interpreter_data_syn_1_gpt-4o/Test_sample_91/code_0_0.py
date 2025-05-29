from itertools import product

# Define the possible outcomes for each die
d20 = range(1, 21)
d10 = range(1, 11)
d4 = range(1, 5)
d3 = range(1, 4)

# Total number of possible outcomes
total_outcomes = 20 * 10 * 4 * 3

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible combinations of dice rolls
for roll in product(d20, d10, d4, d3):
    if sum(roll) >= 21:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)