from fractions import Fraction

# Define the number of sides for each die
sides = [20, 17, 13, 10]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d17 in range(1, 18):
        for d13 in range(1, 14):
            for d10 in range(1, 11):
                if d20 + d17 + d13 + d10 >= 29:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)