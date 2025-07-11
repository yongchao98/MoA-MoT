# Define the two numbers for the problem
a = 6
b = 12

# The discovered pattern is: a * b + a - (b / a)
# We calculate the result using this formula.
# Note that in all given cases, b is a multiple of a, so b / a results in an integer.
result = a * b + a - (b / a)

# The result is a float after division, so we convert it to an integer.
final_result = int(result)

# Print the final equation with each number, as requested.
print(f"{a} + {b} = {final_result}")