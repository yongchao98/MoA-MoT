from sympy import symbols, Eq, solve

# Define the variables
h, h_cable = symbols('h h_cable')

# Equation from similar triangles
eq1 = Eq(h / 3, (1.5 + h_cable) / 2.5)

# Equation for the height of the cable at 2.5 meters
eq2 = Eq(h_cable, h * (2.5 / 3))

# Solve the system of equations
solution = solve((eq1, eq2), (h, h_cable))

# Extract the height of the pole
pole_height = solution[h]
print(pole_height)