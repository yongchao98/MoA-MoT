# A number is a "segmented number" if it cannot be expressed as the sum of
# at least two consecutive positive integers.
# This property holds true only for numbers that are powers of 2.
# The sequence of segmented numbers is 2^0, 2^1, 2^2, 2^3, ...
# The first element is 1 (2^0), the second is 2 (2^1), and so on.
# We need to find the 50th element in this sequence.

# The formula for the n-th element is 2^(n-1).
n = 50

# Calculate the 50th element.
base = 2
exponent = n - 1
result = base ** exponent

# The problem asks to output each number in the final equation.
# The final equation is 2^49 = result.
# The numbers are the base (2), the exponent (49), and the result.
print(f"{base} ** {exponent} = {result}")