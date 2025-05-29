# Define the number of sides for each die
sides = [20, 19, 16, 6]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d19 in range(1, 20):
        for d16 in range(1, 17):
            for d6 in range(1, 7):
                if d20 + d19 + d16 + d6 >= 32:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)