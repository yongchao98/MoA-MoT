from fractions import Fraction

# Define the number of sides for each die
sides = [20, 14, 12, 9]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d14 in range(1, 15):
        for d12 in range(1, 13):
            for d9 in range(1, 10):
                if d20 + d14 + d12 + d9 >= 26:
                    favorable_outcomes += 1

# Calculate the probability as a reduced fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the probability
print(probability)