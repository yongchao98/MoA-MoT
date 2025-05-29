# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d9 in range(1, 10):
        for d7_1 in range(1, 8):
            for d7_2 in range(1, 8):
                # Calculate the sum of the current combination
                total = d20 + d9 + d7_1 + d7_2
                # Check if the sum is 20 or higher
                if total >= 20:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 9 * 7 * 7

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the reduced fraction
print(probability)