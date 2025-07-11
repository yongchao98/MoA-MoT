import math

# After solving the differential equation and simplifying the expression,
# we are left with the final equation to calculate.
# The simplified expression is (3/2) * 10**(10/3) + 37/4.
# The numbers in this final equation are 3, 2, 10, 10, 3, 37, 4.

# Assigning each number to a variable
num1 = 3
den1 = 2
base = 10
pow_num = 10
pow_den = 3
num2 = 37
den2 = 4

# Calculate the result
term1 = (num1 / den1) * (base ** (pow_num / pow_den))
term2 = num2 / den2
result = term1 + term2

# Print the final equation and its value
print(f"The final simplified equation is: ({num1}/{den1}) * {base}**({pow_num}/{pow_den}) + ({num2}/{den2})")
print(f"The value is: {result}")