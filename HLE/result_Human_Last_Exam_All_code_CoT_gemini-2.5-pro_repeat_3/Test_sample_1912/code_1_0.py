# Define the two numbers for the final equation
a = 6
b = 12

# The pattern is determined by the ratio of b to a.
# If b / a == 2, the formula is a * b + (a - 2).
# If b / a == 3, the formula is a * b.

# For 6 + 12, the ratio is 12 / 6 = 2.
# So we use the first formula.
result = a * b + (a - 2)

# Print the final equation with the calculated result.
print(f"{a} + {b} = {result}")