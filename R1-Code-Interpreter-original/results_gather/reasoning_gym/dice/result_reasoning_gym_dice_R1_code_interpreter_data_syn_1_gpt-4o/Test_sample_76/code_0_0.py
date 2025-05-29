from itertools import product

# Define the range of each die
d20 = range(1, 21)
d14 = range(1, 15)
d6 = range(1, 7)
d4 = range(1, 5)

# Total number of possible outcomes
total_outcomes = len(d20) * len(d14) * len(d6) * len(d4)

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible combinations of dice rolls
for roll in product(d20, d14, d6, d4):
    if sum(roll) >= 28:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)