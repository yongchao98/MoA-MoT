# Calculate the probability of drawing two red marbles
probability_first_red = 2 / 15
probability_second_red = 1 / 14

# Combined probability
combined_probability = probability_first_red * probability_second_red

# Print the result as a fraction
from fractions import Fraction
result = Fraction(combined_probability).limit_denominator()
print(result)