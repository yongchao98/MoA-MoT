# Define the numbers for the final equation
a = 6
b = 12

# The pattern discovered is:
# result = a * b + (a - 2) if b is twice a
# result = a * b            otherwise

# Apply the pattern to the given numbers a=6, b=12
# First, check the condition if b is twice a
if b == 2 * a:
    # Calculate the result using the first rule
    result = a * b + (a - 2)
else:
    # Calculate the result using the second rule
    result = a * b

# Print the final equation with the calculated result
# The output format is a + b = result, as seen in the problem description.
print(f"{a} + {b} = {result}")