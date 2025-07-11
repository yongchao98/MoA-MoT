# Define the input numbers for the final equation
a = 6
b = 12

# The discovered pattern is that 'a + b' translates to the formula: a * b + a - (b / a)
# Let's apply this formula to calculate the result.
# Note that integer division '//' can be used since b is always a multiple of a in the sequence.
result = a * b + a - (b // a)

# Print the final equation with the calculated result
print(f"{a} + {b} = {result}")