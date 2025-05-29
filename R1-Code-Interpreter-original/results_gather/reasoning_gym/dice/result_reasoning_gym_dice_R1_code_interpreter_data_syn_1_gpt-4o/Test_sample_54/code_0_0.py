from itertools import product

# Define the number of sides for each die
sides = [20, 19, 12, 2]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Generate all possible outcomes
all_rolls = product(range(1, 21), range(1, 20), range(1, 13), range(1, 3))

# Count the number of favorable outcomes
favorable_outcomes = sum(1 for roll in all_rolls if sum(roll) >= 20)

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes).limit_denominator()

print(probability)