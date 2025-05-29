from sympy import symbols, cos, diff

# Define the variable and the function
x = symbols('x')
integral_result = -1/2 * cos(2*x)

# Differentiate the result
differentiated_result = diff(integral_result, x)

# Output the differentiated result
print(differentiated_result)