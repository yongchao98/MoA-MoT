# Define the variables for the problem.
a = 6
b = 12

# The pattern discovered from the sequence is: a * b + (5*a - b - 6) / 3
# Let's calculate the result based on this pattern.
# The division here must be float division.
result = a * b + (5 * a - b - 6) / 3

# Print the final equation with the calculated result.
# The result is a float, so we convert it to an integer for the final output.
print(f"{a} + {b} = {int(result)}")