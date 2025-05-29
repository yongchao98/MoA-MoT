# Define the number of sides for each die
sides_d20 = 20
sides_d11 = 11
sides_d7 = 7
sides_d5 = 5

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d11 in range(1, sides_d11 + 1):
        for d7 in range(1, sides_d7 + 1):
            for d5 in range(1, sides_d5 + 1):
                if d20 + d11 + d7 + d5 >= 28:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d11 * sides_d7 * sides_d5

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the probability
print(probability)