# Define the two numbers for the final equation
a = 6
b = 12

# The pattern depends on the ratio of b to a.
# If b is twice a, the pattern is a * b + a - 2.
# If b is three times a, the pattern is a * b.
if b / a == 2:
    result = a * b + a - 2
elif b / a == 3:
    result = a * b
else:
    result = "Undefined pattern"

# Print the final equation including the two input numbers and the calculated result.
print(f"{a} + {b} = {result}")