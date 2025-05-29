# Define the range of each die
d20 = range(1, 21)
d9 = range(1, 10)
d7 = range(1, 8)
d4 = range(1, 5)

# Total number of possible outcomes
total_outcomes = 20 * 9 * 7 * 4

# Count the number of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll_d20 in d20:
    for roll_d9 in d9:
        for roll_d7 in d7:
            for roll_d4 in d4:
                if roll_d20 + roll_d9 + roll_d7 + roll_d4 >= 23:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Output the probability as a reduced fraction
print(probability)