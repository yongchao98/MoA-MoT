# Initialize the count of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d5 in range(1, 6):
        for d3 in range(1, 4):
            for d2 in range(1, 3):
                # Calculate the sum of the current combination
                total = d20 + d5 + d3 + d2
                # Check if the sum is 12 or higher
                if total >= 12:
                    successful_outcomes += 1

# Total number of possible outcomes
total_outcomes = 20 * 5 * 3 * 2

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)