# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d15_1 in range(1, 16):
        for d15_2 in range(1, 16):
            for d5 in range(1, 6):
                # Calculate the sum of the current roll
                total = d20 + d15_1 + d15_2 + d5
                # Check if the sum is 23 or higher
                if total >= 23:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = 20 * 15 * 15 * 5

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes).limit_denominator()

# Print the probability
print(probability)