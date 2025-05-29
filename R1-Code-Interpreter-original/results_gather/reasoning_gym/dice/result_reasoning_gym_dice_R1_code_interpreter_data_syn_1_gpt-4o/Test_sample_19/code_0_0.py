from itertools import product

# Define the number of sides for each die
sides = [20, 19, 18, 9]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Count the number of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for roll in product(range(1, sides[0] + 1), range(1, sides[1] + 1), range(1, sides[2] + 1), range(1, sides[3] + 1)):
    if sum(roll) >= 29:
        successful_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)