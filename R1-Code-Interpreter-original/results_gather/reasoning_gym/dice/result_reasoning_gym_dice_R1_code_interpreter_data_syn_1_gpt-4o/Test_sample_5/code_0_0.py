# Define the number of sides for each die
sides_d20 = 20
sides_d18 = 18
sides_d4 = 4
sides_d3 = 3

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d18 in range(1, sides_d18 + 1):
        for d4 in range(1, sides_d4 + 1):
            for d3 in range(1, sides_d3 + 1):
                # Calculate the sum of the dice
                total = d20 + d18 + d4 + d3
                # Check if the sum is 27 or higher
                if total >= 27:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d18 * sides_d4 * sides_d3

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the probability as a reduced fraction
print(probability)