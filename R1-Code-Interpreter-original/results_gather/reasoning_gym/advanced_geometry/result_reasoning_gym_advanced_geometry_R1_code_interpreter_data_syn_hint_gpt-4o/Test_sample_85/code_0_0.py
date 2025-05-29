from sympy import symbols, Eq, solve

# Define the variable
x = symbols('x')

# Equation for the intersection of the altitudes
equation = Eq(2 + 5, (1/5)*(x + 5))

# Solve for x
x_value = solve(equation, x)[0]

# y is already known from the altitude from B
y_value = 2

# Print the orthocenter coordinates
print((round(x_value, 3), round(y_value, 3)))