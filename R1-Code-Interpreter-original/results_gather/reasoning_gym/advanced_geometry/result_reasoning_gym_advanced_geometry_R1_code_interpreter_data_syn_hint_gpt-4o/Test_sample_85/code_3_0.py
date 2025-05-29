from sympy import symbols, Eq, solve

# Define the variable
x = symbols('x')

# Equation for the intersection of the altitudes
equation = Eq((1/5)*x + 6, (-1/7)*x + 34/7)

# Solve for x
x_value = solve(equation, x)[0]

# Substitute x back into one of the altitude equations to find y
y_value = (1/5)*x_value + 6

# Print the orthocenter coordinates
print((round(x_value, 3), round(y_value, 3)))