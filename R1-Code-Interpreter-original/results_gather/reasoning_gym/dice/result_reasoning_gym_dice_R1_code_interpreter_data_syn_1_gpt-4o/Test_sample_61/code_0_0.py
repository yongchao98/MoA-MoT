from fractions import Fraction

# Define the number of sides for each die
sides = [20, 17, 13, 5]

# Calculate the total number of outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d17 in range(1, 18):
        for d13 in range(1, 14):
            for d5 in range(1, 6):
                if d20 + d17 + d13 + d5 >= 23:
                    favorable_outcomes += 1

# Calculate the probability as a reduced fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability
print(probability)