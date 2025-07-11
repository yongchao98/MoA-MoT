import math

# The problem simplifies to calculating the value of the expression:
# (3/2) * 10^(10/3) + 37/4
# We define the numbers in this final equation.
num1 = 3
den1 = 2
base = 10
pow_num = 10
pow_den = 3
num2 = 37
den2 = 4

# Calculate the two terms of the sum
term1 = (num1 / den1) * (base ** (pow_num / pow_den))
term2 = num2 / den2

# Calculate the final result
result = term1 + term2

# Print the final equation with all its numbers and the computed result
print(f"The final calculation is:")
print(f"({num1}/{den1}) * {base}^({pow_num}/{pow_den}) + {num2}/{den2} = {result}")
