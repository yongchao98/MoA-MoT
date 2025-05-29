# Define the number of sides for each die
sides_d20 = 20
sides_d8 = 8
sides_d6 = 6
sides_d4 = 4

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d8 * sides_d6 * sides_d4

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d8 in range(1, sides_d8 + 1):
        for d6 in range(1, sides_d6 + 1):
            for d4 in range(1, sides_d4 + 1):
                if d20 + d8 + d6 + d4 >= 19:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)