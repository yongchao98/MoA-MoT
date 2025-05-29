# Define the range of each die
d20_range = range(1, 21)
d16_range = range(1, 17)
d6_range = range(1, 7)

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in d20_range:
    for d16_1 in d16_range:
        for d16_2 in d16_range:
            for d6 in d6_range:
                if d20 + d16_1 + d16_2 + d6 >= 27:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = 20 * 16 * 16 * 6

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)