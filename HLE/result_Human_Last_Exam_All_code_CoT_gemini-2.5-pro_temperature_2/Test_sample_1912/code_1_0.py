# Assign the values for the final equation
a = 6
b = 12

# The pattern is a * b + a - (b / a)
# In all examples, b is divisible by a, so we can use integer division //.
result = a * b + a - (b // a)

# Print the final equation with the calculated result
print(f"{a} + {b} = {result}")