# The problem asks for the 50th segmented number.
# A number is segmented if it cannot be written as the sum of two or more
# consecutive positive integers.
# A positive integer can be written as a sum of two or more consecutive
# positive integers if and only if it is not a power of 2.
# Therefore, segmented numbers are powers of 2.
# The sequence starts with the first element being 1 (2^0) and the second being 2 (2^1).
# The n-th segmented number is 2^(n-1).

# We need to compute the 50th element of this sequence.
n = 50

# The formula for the n-th element is 2^(n-1).
base = 2
exponent = n - 1

# Calculate the result. Python handles large integers automatically.
result = base ** exponent

# Print the final equation with all its numbers as requested.
print(f"{base} ** {exponent} = {result}")