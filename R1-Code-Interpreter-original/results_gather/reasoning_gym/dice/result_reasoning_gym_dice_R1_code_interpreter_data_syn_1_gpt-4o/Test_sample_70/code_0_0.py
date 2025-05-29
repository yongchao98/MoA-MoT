# Define the number of sides for each die
sides_d20 = 20
sides_d13 = 13
sides_d7 = 7
sides_d6 = 6

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d13 in range(1, sides_d13 + 1):
        for d7 in range(1, sides_d7 + 1):
            for d6 in range(1, sides_d6 + 1):
                # Calculate the sum of the dice
                total = d20 + d13 + d7 + d6
                # Check if the sum is 26 or higher
                if total >= 26:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d13 * sides_d7 * sides_d6

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)