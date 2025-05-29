# Simulate rolling two dice and calculate the expected number of 1's
from itertools import product

# All possible outcomes for two dice
outcomes = list(product(range(1, 7), repeat=2))

# Calculate the expected number of 1's
expected_ones = sum((die1 == 1) + (die2 == 1) for die1, die2 in outcomes) / len(outcomes)

print(expected_ones)