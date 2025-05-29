# Initialize the count of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d7_1 in range(1, 8):
        for d7_2 in range(1, 8):
            for d4 in range(1, 5):
                # Calculate the sum of the current combination
                total = d20 + d7_1 + d7_2 + d4
                # Check if the sum is 21 or higher
                if total >= 21:
                    successful_outcomes += 1

# Total number of possible outcomes
total_outcomes = 20 * 7 * 7 * 4

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)