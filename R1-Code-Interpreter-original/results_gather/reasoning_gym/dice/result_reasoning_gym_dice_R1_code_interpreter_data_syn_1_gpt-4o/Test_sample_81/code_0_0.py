from fractions import Fraction

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d10 in range(1, 11):
        for d2_1 in range(1, 3):
            for d2_2 in range(1, 3):
                total = d20 + d10 + d2_1 + d2_2
                if total >= 21:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 10 * 2 * 2

# Calculate the probability as a fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)