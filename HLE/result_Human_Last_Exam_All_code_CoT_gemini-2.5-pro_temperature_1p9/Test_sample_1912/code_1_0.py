# Define the input numbers for the final equation
a = 6
b = 12

# Apply the discovered pattern: a * b + (a - b / a)
# Using integer division // as all intermediate results are integers
result = a * b + (a - b // a)

# Print the final equation with each number
print(f"{a} + {b} = {result}")