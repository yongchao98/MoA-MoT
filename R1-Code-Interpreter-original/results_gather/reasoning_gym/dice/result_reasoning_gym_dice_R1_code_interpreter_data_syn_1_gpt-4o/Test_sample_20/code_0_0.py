# Define the range of each die
d20 = range(1, 21)
d18 = range(1, 19)
d10 = range(1, 11)
d3 = range(1, 4)

# Initialize counters
total_outcomes = 0
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll_d20 in d20:
    for roll_d18 in d18:
        for roll_d10 in d10:
            for roll_d3 in d3:
                total_outcomes += 1
                if roll_d20 + roll_d18 + roll_d10 + roll_d3 >= 34:
                    favorable_outcomes += 1

# Calculate the probability as a reduced fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the result
print(probability)