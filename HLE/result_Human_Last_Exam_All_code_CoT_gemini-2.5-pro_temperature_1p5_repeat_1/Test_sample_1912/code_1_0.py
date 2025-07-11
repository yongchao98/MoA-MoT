# Define the input numbers for the equation we need to solve
a = 6
b = 12

# The pattern discovered is: a * b + a - (b / a)
# Let's apply this pattern to solve for '?'
multiplication_part = a * b
division_part = b // a # Using integer division as per the pattern
addition_part = a
subtraction_part = division_part

# Calculate the final result
result = multiplication_part + addition_part - subtraction_part

# Print the final equation as requested, showing each number
print(f"{a} + {b} = {result}")