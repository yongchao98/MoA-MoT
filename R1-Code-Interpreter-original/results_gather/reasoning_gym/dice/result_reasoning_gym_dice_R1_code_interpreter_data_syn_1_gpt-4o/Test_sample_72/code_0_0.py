# Initialize the count of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d9 in range(1, 10):
        for d8 in range(1, 9):
            for d2 in range(1, 3):
                # Calculate the sum of the current roll
                total = d20 + d9 + d8 + d2
                # Check if the sum is 21 or higher
                if total >= 21:
                    successful_outcomes += 1

# Total number of possible outcomes
total_outcomes = 20 * 9 * 8 * 2

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)