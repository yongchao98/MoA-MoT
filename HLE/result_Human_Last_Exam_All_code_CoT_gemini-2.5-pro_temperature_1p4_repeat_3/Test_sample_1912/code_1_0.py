# Define the numbers for the final equation
a = 6
b = 12

# Apply the discovered pattern: result = a * b + a - (b / a)
# We use integer division // as the relationship is always an integer in the examples
result = a * b + a - (b // a)

# Print the final equation with the calculated result
print(f"{a} + {b} = {result}")