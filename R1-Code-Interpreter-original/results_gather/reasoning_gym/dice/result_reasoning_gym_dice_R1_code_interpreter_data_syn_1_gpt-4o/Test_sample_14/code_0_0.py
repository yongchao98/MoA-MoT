# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d14 in range(1, 15):
        for d4 in range(1, 5):
            for d3 in range(1, 4):
                if d20 + d14 + d4 + d3 >= 25:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 14 * 4 * 3

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)