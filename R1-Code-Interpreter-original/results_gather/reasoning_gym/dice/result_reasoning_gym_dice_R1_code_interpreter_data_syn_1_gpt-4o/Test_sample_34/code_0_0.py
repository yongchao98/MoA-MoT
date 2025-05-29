from itertools import product

# Define the number of sides for each die
sides = [20, 15, 11, 4]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, sides[0] + 1), range(1, sides[1] + 1), range(1, sides[2] + 1), range(1, sides[3] + 1)):
    if sum(roll) >= 26:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)