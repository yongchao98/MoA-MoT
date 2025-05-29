from itertools import product

# Define the number of sides for each die
sides = [20, 16, 9, 4]

# Calculate total number of outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Generate all possible outcomes
all_outcomes = product(range(1, sides[0] + 1), range(1, sides[1] + 1), range(1, sides[2] + 1), range(1, sides[3] + 1))

# Count the number of favorable outcomes
favorable_outcomes = sum(1 for outcome in all_outcomes if sum(outcome) >= 19)

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)