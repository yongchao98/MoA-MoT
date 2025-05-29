# Initialize the count of successful outcomes
successful_outcomes = 0

# Iterate over all possible outcomes
for d20 in range(1, 21):
    for d18 in range(1, 19):
        for d10 in range(1, 11):
            for d5 in range(1, 6):
                # Calculate the sum of the current roll
                total = d20 + d18 + d10 + d5
                # Check if the sum is 35 or higher
                if total >= 35:
                    successful_outcomes += 1

# Total number of possible outcomes
total_outcomes = 20 * 18 * 10 * 5

# Calculate the probability as a fraction
from fractions import Fraction
probability = Fraction(successful_outcomes, total_outcomes)

# Print the reduced fraction
print(probability)