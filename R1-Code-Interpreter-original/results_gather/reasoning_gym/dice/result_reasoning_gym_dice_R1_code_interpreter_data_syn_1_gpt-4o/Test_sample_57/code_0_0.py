# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d10 in range(1, 11):
        for d6 in range(1, 7):
            for d4 in range(1, 5):
                # Calculate the sum of the current combination
                total = d20 + d10 + d6 + d4
                # Check if the sum is 21 or higher
                if total >= 21:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, 4800)

# Print the reduced fraction
print(probability)