from fractions import Fraction

# Total number of outcomes
total_outcomes = 20 * 13 * 5 * 2

# Count favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d13 in range(1, 14):
        for d5 in range(1, 6):
            for d2 in range(1, 3):
                if d20 + d13 + d5 + d2 >= 20:
                    favorable_outcomes += 1

# Calculate probability as a reduced fraction
probability = Fraction(favorable_outcomes, total_outcomes)

print(probability)