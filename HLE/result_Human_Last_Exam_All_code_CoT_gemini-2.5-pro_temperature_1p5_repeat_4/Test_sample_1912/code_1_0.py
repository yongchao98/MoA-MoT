# Define the input numbers for the final equation.
a = 6
b = 12

# Apply the discovered pattern: result = a * b + a - (b / a)
# Division (b / a) should be performed first.
# For integer results, we can use integer division //.
result = a * b + a - (b / a)

# Print the final equation with the calculated result.
# The result is a float from the division, so we cast it to an integer.
print(f"{a} + {b} = {int(result)}")