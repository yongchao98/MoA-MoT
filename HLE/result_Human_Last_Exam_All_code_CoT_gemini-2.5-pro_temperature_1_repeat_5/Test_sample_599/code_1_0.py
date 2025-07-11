# A segmented number is a positive integer that cannot be written
# as the sum of two or more consecutive positive integers.
# Mathematical analysis shows that this property holds true only for numbers that are powers of 2.
# The sequence of segmented numbers is 1, 2, 4, 8, 16, ... which corresponds to 2^0, 2^1, 2^2, 2^3, 2^4, ...
# The n-th term in this sequence follows the formula: 2^(n-1).

# We are asked to find the 50th term in this sequence.
n = 50

# The final equation to compute the result is: result = 2^(n-1).
# For n=50, the equation becomes result = 2^49.
# The numbers in this final equation are the base (2), the exponent (49), and the result itself.
base = 2
exponent = n - 1

# Calculate the 50th segmented number
result = base ** exponent

# Print the final answer.
print(result)