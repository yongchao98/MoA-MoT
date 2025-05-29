total_outcomes = 0
successful_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d11 in range(1, 12):
        for d9_1 in range(1, 10):
            for d9_2 in range(1, 10):
                total_outcomes += 1
                if d20 + d11 + d9_1 + d9_2 >= 24:
                    successful_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

print(probability)