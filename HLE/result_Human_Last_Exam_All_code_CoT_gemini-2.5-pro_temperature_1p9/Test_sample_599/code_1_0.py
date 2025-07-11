# A segmented number is a positive integer that cannot be expressed as the sum of two or more
# consecutive positive integers. It can be shown that these numbers are precisely the powers of 2.
# The sequence is 2^0, 2^1, 2^2, 2^3, ...
# We need to find the 50th element of this sequence.
# The n-th element is given by the formula 2^(n-1).
# For the 50th element, n = 50.
# So, we need to calculate 2^(50-1), which is 2^49.

# The exponent in the equation
n = 50
exponent = n - 1

# The base of the power
base = 2

# Calculate the result
result = base ** exponent

# The final equation includes the base (2), the exponent (49), and the result.
# We print all these numbers as requested.
print(f"{base} ** {exponent} = {result}")