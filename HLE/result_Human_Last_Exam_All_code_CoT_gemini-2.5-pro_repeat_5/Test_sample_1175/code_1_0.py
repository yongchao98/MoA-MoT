import math

# Define the numbers in the final simplified equation
num1 = 3
den1 = 2
base = 10
exp_num = 10
exp_den = 3
num2 = 37
den2 = 4

# Calculate the final result
term1 = (num1 / den1) * (base ** (exp_num / exp_den))
term2 = num2 / den2
result = term1 + term2

# Print the final equation using the defined numbers
print(f"The final simplified expression to calculate is:")
print(f"({num1}/{den1}) * {base}**({exp_num}/{exp_den}) + ({num2}/{den2})")
print(f"The calculated value is: {result}")