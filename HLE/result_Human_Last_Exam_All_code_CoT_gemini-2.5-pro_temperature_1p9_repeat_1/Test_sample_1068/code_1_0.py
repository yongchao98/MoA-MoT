import math

# The value of n as per the problem description.
n = 13
# The value of p for which we need to calculate l(p).
p = 14

# After simplifying the complex definitions, the function l(p) can be
# expressed as l(p) = p * log(4).
# For p = 14, the value is 14 * log(4).
# This can also be written as 14 * 2 * log(2) = 28 * log(2).

# Calculate the result
result = 28 * math.log(2)

# Print the final equation with all numbers, as requested.
# The formula is l(14) = 28 * log(2)
print("l(14) =", 28, "*", "log(2)", "=", result)