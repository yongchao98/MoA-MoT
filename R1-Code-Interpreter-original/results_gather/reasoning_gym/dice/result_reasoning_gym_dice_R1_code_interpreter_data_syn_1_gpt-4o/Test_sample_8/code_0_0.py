# Define the number of sides for each die
sides_d20 = 20
sides_d18 = 18
sides_d16 = 16
sides_d3 = 3

# Initialize the count of favorable outcomes
favorable_outcomes = 0

# Iterate over all possible outcomes
for roll_d20 in range(1, sides_d20 + 1):
    for roll_d18 in range(1, sides_d18 + 1):
        for roll_d16 in range(1, sides_d16 + 1):
            for roll_d3 in range(1, sides_d3 + 1):
                # Calculate the sum of the current roll
                total = roll_d20 + roll_d18 + roll_d16 + roll_d3
                # Check if the total is 28 or higher
                if total >= 28:
                    favorable_outcomes += 1

# Calculate the total number of possible outcomes
total_outcomes = sides_d20 * sides_d18 * sides_d16 * sides_d3

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(favorable_outcomes, total_outcomes)

# Print the probability as a reduced fraction
print(probability)