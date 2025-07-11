# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be written as the sum
# of two or more consecutive smaller positive integers.
# It is a mathematical fact that these numbers are the powers of 2.
# The sequence of segmented numbers is: 2^0, 2^1, 2^2, 2^3, ...
# So, the n-th segmented number is 2^(n-1).

# We need to find the 50th element in this sequence.
n = 50
base = 2
exponent = n - 1

# Calculate the result
result = base ** exponent

# Print the final equation with all its numbers
print(f"{base} ** {exponent} = {result}")