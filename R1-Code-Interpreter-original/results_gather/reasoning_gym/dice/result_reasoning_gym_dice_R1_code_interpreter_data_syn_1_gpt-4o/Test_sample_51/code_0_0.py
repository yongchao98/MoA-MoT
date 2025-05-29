from itertools import product

# Define the number of sides for each die
sides = [20, 17, 16, 4]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for outcome in product(range(1, 21), range(1, 18), range(1, 17), range(1, 5)):
    if sum(outcome) >= 25:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability
print(probability)