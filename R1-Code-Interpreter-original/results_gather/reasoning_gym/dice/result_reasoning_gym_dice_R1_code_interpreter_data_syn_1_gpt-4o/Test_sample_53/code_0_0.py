# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d12 in range(1, 13):
        for d8 in range(1, 9):
            for d5 in range(1, 6):
                if d20 + d12 + d8 + d5 >= 31:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 12 * 8 * 5

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)