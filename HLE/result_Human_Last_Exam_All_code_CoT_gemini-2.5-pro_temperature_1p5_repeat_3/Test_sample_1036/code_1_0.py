# The problem is to find the number of positive integers n <= lcm(1, ..., 100)
# such that n gives distinct remainders when divided by k = 2, 3, ..., 100.
# As derived from the analysis, the number of ways to choose a valid sequence
# of remainders (r_2, r_3, ..., r_100) that satisfy the given conditions is 2^99.
#
# Each such sequence corresponds to a unique integer n in the specified range
# due to the Chinese Remainder Theorem.
# Therefore, we need to compute the value of 2^99.

# The base of the power
base = 2

# The exponent
exponent = 99

# Calculate the result
result = base**exponent

# Print the final equation and the result
print("The number of such integers is given by the equation:")
print(base, "^", exponent, "=", result)