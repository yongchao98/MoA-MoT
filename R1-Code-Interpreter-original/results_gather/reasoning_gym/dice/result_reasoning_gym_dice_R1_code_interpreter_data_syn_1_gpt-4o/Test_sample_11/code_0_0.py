# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d11 in range(1, 12):
        for d6 in range(1, 7):
            for d3 in range(1, 4):
                if d20 + d11 + d6 + d3 >= 23:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 11 * 6 * 3

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)