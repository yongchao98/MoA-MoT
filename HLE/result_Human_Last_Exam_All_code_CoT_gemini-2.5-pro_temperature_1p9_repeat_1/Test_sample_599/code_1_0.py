# The problem asks for the 50th segmented number.
# A segmented number is a positive integer that cannot be expressed
# as the sum of two or more consecutive positive integers.
# Through mathematical analysis, it can be shown that these numbers
# are precisely the powers of 2.
# The sequence starts with 2^0=1, 2^1=2, 2^2=4, and so on.

# The nth segmented number is given by the formula 2^(n-1).
# For the 50th element, n=50.

# The position in the sequence
n = 50

# The exponent is n-1
exponent = n - 1

# The base is 2
base = 2

# Calculate the result
result = base ** exponent

# Print the final equation with all numbers
print(f"The 50th segmented number is the result of the equation:")
print(f"{base} ** {exponent} = {result}")