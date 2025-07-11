# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed as the sum
# of two or more consecutive smaller positive integers.

# Through number theory, it can be proven that these numbers are precisely the powers of 2.
# A number N can be written as a sum of k>=2 consecutive integers if and only if N is not a power of 2.
# Therefore, the numbers that CANNOT be written this way are the powers of 2.

# The sequence of segmented numbers is:
# 1st element: 2^0 = 1
# 2nd element: 2^1 = 2
# 3rd element: 2^2 = 4
# ...
# n-th element: 2^(n-1)

# We need to find the 50th element of this sequence.
n = 50

# The formula for the n-th element is 2^(n-1).
base = 2
exponent = n - 1

# Calculate the result of the equation.
result = base ** exponent

# As requested, we output the numbers in the final equation, which are
# the base and the exponent, followed by the result.
print(f"The 50th segmented number is the result of the equation: {base} ** {exponent}")
print(f"Result: {result}")