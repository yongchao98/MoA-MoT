from itertools import product

# Define the number of sides for each die
sides = [20, 11, 8, 5]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, 21), range(1, 12), range(1, 9), range(1, 6)):
    if sum(roll) >= 25:
        favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)