# Define the variables for the problem
a = 6
b = 12

# The pattern is a * b + a - (b / a)
# In Python, we use // for integer division
result = a * b + a - (b // a)

# Print the final equation with the calculated result, showing the steps of the pattern.
print(f"{a} + {b} = {a} * {b} + ({a} - {b} // {a}) = {result}")
