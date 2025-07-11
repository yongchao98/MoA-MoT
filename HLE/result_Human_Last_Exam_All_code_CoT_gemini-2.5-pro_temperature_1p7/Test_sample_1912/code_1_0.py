# Define the variables for the final equation in the sequence.
a = 6
b = 12

# The pattern discovered from the sequence is: result = a * b + a - (b / a).
# We apply this formula to find the missing number.
# We use integer division // since all previous results are integers.
result = a * b + a - (b // a)

# Print the complete final equation, including the numbers a, b, and the calculated result.
print(f"{a} + {b} = {result}")