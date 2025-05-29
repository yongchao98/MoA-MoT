# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes for each die
for d20 in range(1, 21):
    for d8 in range(1, 9):
        for d6 in range(1, 7):
            for d4 in range(1, 5):
                # Calculate the sum of the dice
                total = d20 + d8 + d6 + d4
                # Check if the sum is 26 or higher
                if total >= 26:
                    favorable_outcomes += 1

# Total number of possible outcomes
total_outcomes = 20 * 8 * 6 * 4

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)