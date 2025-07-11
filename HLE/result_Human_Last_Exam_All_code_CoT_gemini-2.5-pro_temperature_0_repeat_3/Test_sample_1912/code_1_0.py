# Define the numbers for the final equation
a = 6
b = 12

# The pattern is determined by the ratio of the two numbers.
# If b / a == 3, the result is a * b.
# If b / a == 2, the result is a * b + a - 2.

# For the given numbers, a=6 and b=12, the ratio is 12 / 6 = 2.
# So we apply the second rule.
result = a * b + a - 2

# Print the final equation including the numbers and the result.
print(f"{a} + {b} = {result}")