# The problem is to find the 50th segmented number.
# A segmented number is a positive integer that cannot be written as the
# sum of two or more consecutive smaller positive integers.
# It can be shown that a number has this property if and only if it is a power of 2.
# The sequence of segmented numbers is: 2^0, 2^1, 2^2, 2^3, ...
# The Nth number in this sequence is 2^(N-1).
# We need to find the 50th number.

# The position in the sequence
n = 50

# The base of the power
base = 2

# The exponent for the Nth term is N-1
exponent = n - 1

# Calculate the 50th segmented number
result = base ** exponent

# Print the final equation and the result.
print(f"The {n}th segmented number is the result of the equation: {base} ** {exponent}")
print(f"The result is: {result}")