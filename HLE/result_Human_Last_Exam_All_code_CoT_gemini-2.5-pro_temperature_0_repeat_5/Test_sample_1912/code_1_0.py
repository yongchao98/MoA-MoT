# Define the numbers for the final equation in the sequence
a = 6
b = 12

# Calculate the result using the discovered pattern: a * b + a - (b / a)
# The result must be an integer, as in the examples.
result = int(a * b + a - (b / a))

# The final equation is '6 + 12 = ?'.
# The following print statement outputs each number (6, 12, and the result 76)
# in the format of the final equation.
print(f"{a} + {b} = {result}")