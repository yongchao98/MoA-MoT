# Define the number of sides for each die
sides_d20 = 20
sides_d19 = 19
sides_d17 = 17
sides_d2 = 2

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, sides_d20 + 1):
    for d19 in range(1, sides_d19 + 1):
        for d17 in range(1, sides_d17 + 1):
            for d2 in range(1, sides_d2 + 1):
                # Calculate the sum of the current outcome
                total = d20 + d19 + d17 + d2
                # Check if the sum is 38 or higher
                if total >= 38:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d19 * sides_d17 * sides_d2

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)