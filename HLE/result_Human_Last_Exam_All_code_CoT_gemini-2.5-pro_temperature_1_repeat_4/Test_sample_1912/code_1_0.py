# Define the input numbers for the final equation in the sequence.
a = 6
b = 12

# The pattern discovered from the sequence is: a * b + a - (b / a)
# Here, we apply this pattern to the given numbers.

# Calculate the result using the formula.
# The result needs to be an integer as seen in the examples.
result = int(a * b + a - (b / a))

# Print the full calculation to show how each number is used,
# as requested by the prompt.
print(f"The pattern is: a * b + a - (b / a)")
print(f"For the numbers {a} and {b}, the calculation is:")
print(f"{a} * {b} + {a} - ({b} / {a}) = {a * b} + {a} - {int(b/a)} = {result}")
print(f"So, the final equation is:")
print(f"{a} + {b} = {result}")