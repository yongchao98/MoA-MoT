# Define the number of sides for each die
sides_d20 = 20
sides_d16 = 16
sides_d10 = 10
sides_d9 = 9

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d16 in range(1, sides_d16 + 1):
        for d10 in range(1, sides_d10 + 1):
            for d9 in range(1, sides_d9 + 1):
                if d20 + d16 + d10 + d9 >= 35:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d16 * sides_d10 * sides_d9

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)