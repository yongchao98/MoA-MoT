# Define the number of sides for each die
sides_d20 = 20
sides_d16 = 16
sides_d8 = 8
sides_d4 = 4

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d16 in range(1, sides_d16 + 1):
        for d8 in range(1, sides_d8 + 1):
            for d4 in range(1, sides_d4 + 1):
                # Calculate the sum of the current combination
                total = d20 + d16 + d8 + d4
                # Check if the sum is 21 or higher
                if total >= 21:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d16 * sides_d8 * sides_d4

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)