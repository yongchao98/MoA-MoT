# The problem asks for the 50th segmented number.
# Based on the definition, segmented numbers are those that cannot be written
# as the sum of two or more consecutive positive integers.
# This sequence is equivalent to the powers of 2: 2^0, 2^1, 2^2, etc.
# The n-th term in this sequence is 2^(n-1).

# For the 50th term, n is 50.
n = 50

# The equation is base ^ exponent = result.
base = 2
exponent = n - 1

# Calculate the result.
result = base ** exponent

# Per the instructions, we print all numbers involved in the final equation.
# The final equation is 2^49 = result.
print(f"{base} ** {exponent} = {result}")