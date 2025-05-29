# Define the number of sides for each die
sides_d20 = 20
sides_d12 = 12
sides_d9 = 9
sides_d5 = 5

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d12 in range(1, sides_d12 + 1):
        for d9 in range(1, sides_d9 + 1):
            for d5 in range(1, sides_d5 + 1):
                if d20 + d12 + d9 + d5 >= 28:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d12 * sides_d9 * sides_d5

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the probability as a reduced fraction
print(probability)