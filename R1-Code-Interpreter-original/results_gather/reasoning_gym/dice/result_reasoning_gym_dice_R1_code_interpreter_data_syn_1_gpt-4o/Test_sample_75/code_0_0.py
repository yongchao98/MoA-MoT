from itertools import product

# Define the number of sides for each die
sides = [20, 17, 15, 9]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Calculate the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, 21), range(1, 18), range(1, 16), range(1, 10)):
    if sum(roll) >= 32:
        favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)