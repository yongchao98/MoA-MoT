# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the
# sum of two or more consecutive positive integers.
# A number has this property if and only if it is a power of 2.
# So, the sequence of segmented numbers is: 2^0, 2^1, 2^2, ... which is 1, 2, 4, ...

# The n-th term in this sequence is given by the formula 2^(n-1).
# We need to find the 50th term.
n = 50
base = 2

# The formula for the 50th term is base^(n - 1).
subtrahend = 1
exponent = n - subtrahend
result = base ** exponent

# Print the final equation with each of its numbers, as requested.
print(f"To find the {n}th segmented number, we use the formula: {base}^({n} - {subtrahend})")
print(f"This simplifies to the equation: {base}^{exponent}")
print(f"The result is: {result}")