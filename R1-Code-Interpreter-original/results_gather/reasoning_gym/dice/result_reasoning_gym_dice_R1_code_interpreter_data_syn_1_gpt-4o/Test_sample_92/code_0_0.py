from fractions import Fraction

# Define the number of sides for each die
sides = [20, 19, 10, 2]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d19 in range(1, 20):
        for d10 in range(1, 11):
            for d2 in range(1, 3):
                if d20 + d19 + d10 + d2 >= 24:
                    favorable_outcomes += 1

# Calculate the probability as a reduced fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the probability
print(probability)