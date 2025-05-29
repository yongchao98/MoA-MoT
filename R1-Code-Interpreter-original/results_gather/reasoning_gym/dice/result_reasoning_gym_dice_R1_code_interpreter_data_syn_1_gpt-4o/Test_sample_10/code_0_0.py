# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d10 in range(1, 11):
        for d5 in range(1, 6):
            for d2 in range(1, 3):
                if d20 + d10 + d5 + d2 >= 18:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 10 * 5 * 2

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)