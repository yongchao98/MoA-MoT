from fractions import Fraction

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d7 in range(1, 8):
        for d5 in range(1, 6):
            for d3 in range(1, 4):
                if d20 + d7 + d5 + d3 >= 15:
                    favorable_outcomes += 1

# Total possible outcomes
total_outcomes = 20 * 7 * 5 * 3

# Calculate the probability as a fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)