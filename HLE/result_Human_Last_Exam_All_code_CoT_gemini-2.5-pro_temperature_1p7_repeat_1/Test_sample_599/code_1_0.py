# The problem is to find the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed
# as the sum of two or more consecutive smaller positive integers.
# Through mathematical analysis, it can be shown that a number has this property
# if and only if it is a power of 2 (e.g., 1, 2, 4, 8, 16, ...).
# The sequence of segmented numbers is 2^0, 2^1, 2^2, and so on.

# The nth term in this sequence is given by the formula 2^(n-1).
# We need to find the 50th term, so we set n = 50.
# The 50th term is 2^(50-1), which simplifies to 2^49.

# The following code calculates this value and prints the final equation.

# The base of the power
base = 2
# The exponent for the 50th term is 50 - 1
exponent = 49

# Calculate the result
result = base ** exponent

# Print the final equation with all its components
print(f"{base} ** {exponent} = {result}")
