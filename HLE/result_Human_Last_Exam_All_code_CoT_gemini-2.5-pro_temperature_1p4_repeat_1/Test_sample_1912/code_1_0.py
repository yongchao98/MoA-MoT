# Define the numbers for the final equation in the sequence
a = 6
b = 12

# The pattern discovered from the sequence is: result = a * b + a - (b / a)
# This code applies the pattern to solve for '?' in '6 + 12 = ?'

# Calculate the result based on the pattern
# 6 * 12 + 6 - (12 / 6)
# 72 + 6 - 2
# 76
result = a * b + a - (b / a)

# Print the final equation, showing how the numbers are used in the pattern to get the result.
# The result of the division is a float, so we cast the final result to an integer.
print(f"{a} + {b} = {a} * {b} + {a} - {b} / {a} = {int(result)}")