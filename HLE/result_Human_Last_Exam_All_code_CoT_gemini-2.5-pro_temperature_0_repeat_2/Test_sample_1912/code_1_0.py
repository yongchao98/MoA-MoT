# Define the numbers for the final equation
a = 6
b = 12

# The discovered pattern is: result = a * b + a - (b/a)
# First, calculate the ratio k = b/a
k = b / a

# Then, calculate the final result using the pattern
result = a * b + a - k

# Print the complete final equation with the calculated result
# The result is cast to an integer as the pattern consistently produces whole numbers.
print(f"{a} + {b} = {int(result)}")