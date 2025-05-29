# Total number of sides for each die
sides_d20 = 20
sides_d16 = 16
sides_d13 = 13
sides_d6 = 6

# Total number of possible outcomes
total_outcomes = sides_d20 * sides_d16 * sides_d13 * sides_d6

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d16 in range(1, sides_d16 + 1):
        for d13 in range(1, sides_d13 + 1):
            for d6 in range(1, sides_d6 + 1):
                if d20 + d16 + d13 + d6 >= 26:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability as a reduced fraction
print(probability)