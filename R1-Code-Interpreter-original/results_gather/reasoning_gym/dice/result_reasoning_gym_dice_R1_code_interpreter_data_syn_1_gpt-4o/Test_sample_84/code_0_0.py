# Define the number of sides for each die
sides_d20 = 20
sides_d19 = 19
sides_d11 = 11
sides_d2 = 2

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d19 * sides_d11 * sides_d2

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll_d20 in range(1, sides_d20 + 1):
    for roll_d19 in range(1, sides_d19 + 1):
        for roll_d11 in range(1, sides_d11 + 1):
            for roll_d2 in range(1, sides_d2 + 1):
                # Calculate the sum of the current roll
                total_roll = roll_d20 + roll_d19 + roll_d11 + roll_d2
                # Check if the sum is 23 or higher
                if total_roll >= 23:
                    favorable_outcomes += 1

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction representing the probability
print(probability)