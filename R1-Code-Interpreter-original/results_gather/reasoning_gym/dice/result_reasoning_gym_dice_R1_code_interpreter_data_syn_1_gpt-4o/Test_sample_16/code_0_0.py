from itertools import product

# Define the range of each die
d20 = range(1, 21)
d19 = range(1, 20)
d12 = range(1, 13)
d6 = range(1, 7)

# Total number of possible outcomes
total_outcomes = len(d20) * len(d19) * len(d12) * len(d6)

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible combinations of dice rolls
for roll in product(d20, d19, d12, d6):
    if sum(roll) >= 38:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)