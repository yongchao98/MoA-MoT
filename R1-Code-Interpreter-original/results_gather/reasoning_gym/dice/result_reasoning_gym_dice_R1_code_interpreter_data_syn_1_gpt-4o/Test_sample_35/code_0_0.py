# Define the number of sides for each die
sides_d20 = 20
sides_d15 = 15
sides_d5 = 5
sides_d4 = 4

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d15 in range(1, sides_d15 + 1):
        for d5 in range(1, sides_d5 + 1):
            for d4 in range(1, sides_d4 + 1):
                if d20 + d15 + d5 + d4 >= 24:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d15 * sides_d5 * sides_d4

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)