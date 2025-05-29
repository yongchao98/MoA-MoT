# Define the number of sides for each die
sides_d20 = 20
sides_d19 = 19
sides_d12 = 12
sides_d6 = 6

# Total number of possible outcomes
total_outcomes = sides_d20 * sides_d19 * sides_d12 * sides_d6

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d19 in range(1, sides_d19 + 1):
        for d12 in range(1, sides_d12 + 1):
            for d6 in range(1, sides_d6 + 1):
                if d20 + d19 + d12 + d6 >= 31:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)