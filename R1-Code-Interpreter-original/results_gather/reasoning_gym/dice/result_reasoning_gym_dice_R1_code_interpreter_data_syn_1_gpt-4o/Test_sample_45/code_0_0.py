# Define the number of sides for each die
sides_d20 = 20
sides_d12 = 12
sides_d10 = 10
sides_d8 = 8

# Total number of possible outcomes
total_outcomes = sides_d20 * sides_d12 * sides_d10 * sides_d8

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d12 in range(1, sides_d12 + 1):
        for d10 in range(1, sides_d10 + 1):
            for d8 in range(1, sides_d8 + 1):
                if d20 + d12 + d10 + d8 >= 26:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)