from fractions import Fraction

# Define the number of sides for each die
sides = [20, 18, 16, 12]

# Calculate the total number of possible outcomes
total_outcomes = 1
for side in sides:
    total_outcomes *= side

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d18 in range(1, 19):
        for d16 in range(1, 17):
            for d12 in range(1, 13):
                if d20 + d18 + d16 + d12 >= 41:
                    favorable_outcomes += 1

# Calculate the probability as a reduced fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability
print(probability)