from itertools import product

# Define the number of sides for each die
sides = [20, 17, 8, 4]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Count the number of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, 21), range(1, 18), range(1, 9), range(1, 5)):
    if sum(roll) >= 28:
        successful_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

print(probability)