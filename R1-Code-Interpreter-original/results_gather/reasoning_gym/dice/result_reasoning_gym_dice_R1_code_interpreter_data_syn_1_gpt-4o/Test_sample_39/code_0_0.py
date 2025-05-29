# Define the number of sides for each die
sides_d20 = 20
sides_d18 = 18
sides_d12 = 12
sides_d2 = 2

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll_d20 in range(1, sides_d20 + 1):
    for roll_d18 in range(1, sides_d18 + 1):
        for roll_d12 in range(1, sides_d12 + 1):
            for roll_d2 in range(1, sides_d2 + 1):
                # Calculate the sum of the current roll
                total = roll_d20 + roll_d18 + roll_d12 + roll_d2
                # Check if the total is 28 or higher
                if total >= 28:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d18 * sides_d12 * sides_d2

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)