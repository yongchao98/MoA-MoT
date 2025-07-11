# The segmented numbers are powers of 2 (1, 2, 4, 8, ...).
# The n-th segmented number is given by the formula 2^(n-1).
# To find the 50th element, we need to calculate 2^(50-1), which is 2^49.
# The final equation is result = 2 ** 49. The numbers in the equation are
# the base (2), the exponent (49), and the result.

base = 2
exponent = 49
result = base ** exponent

# Print the final equation with all its numbers.
# This format displays each number from the equation "2 ** 49 = result".
print(f"{base} ** {exponent} = {result}")