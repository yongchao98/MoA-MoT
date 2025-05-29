# Initialize variables for the number of sides on each die
sides_d20 = 20
sides_d19_1 = 19
sides_d19_2 = 19
sides_d18 = 18

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll_d20 in range(1, sides_d20 + 1):
    for roll_d19_1 in range(1, sides_d19_1 + 1):
        for roll_d19_2 in range(1, sides_d19_2 + 1):
            for roll_d18 in range(1, sides_d18 + 1):
                # Calculate the sum of the current roll
                total = roll_d20 + roll_d19_1 + roll_d19_2 + roll_d18
                # Check if the sum is 32 or higher
                if total >= 32:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d19_1 * sides_d19_2 * sides_d18

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)