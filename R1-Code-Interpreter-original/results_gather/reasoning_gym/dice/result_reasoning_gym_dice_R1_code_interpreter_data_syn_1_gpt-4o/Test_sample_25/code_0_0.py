# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d13 in range(1, 14):
        for d3 in range(1, 4):
            for d2 in range(1, 3):
                # Calculate the sum of the current combination
                total = d20 + d13 + d3 + d2
                # Check if the sum is 18 or higher
                if total >= 18:
                    favorable_outcomes += 1

# Total number of possible outcomes
total_outcomes = 20 * 13 * 3 * 2

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)