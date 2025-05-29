# Define the number of sides for each die
sides_d20 = 20
sides_d11 = 11
sides_d8 = 8
sides_d5 = 5

# Total number of possible outcomes
total_outcomes = sides_d20 * sides_d11 * sides_d8 * sides_d5

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d11 in range(1, sides_d11 + 1):
        for d8 in range(1, sides_d8 + 1):
            for d5 in range(1, sides_d5 + 1):
                if d20 + d11 + d8 + d5 >= 27:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability as a reduced fraction
print(probability)