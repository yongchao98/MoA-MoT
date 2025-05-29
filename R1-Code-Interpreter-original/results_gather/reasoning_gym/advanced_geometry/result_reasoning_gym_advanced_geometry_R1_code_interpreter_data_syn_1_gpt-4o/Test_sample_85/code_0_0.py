from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Equation of the altitude from C
altitude_C = Eq(y, -1/7 * x + 54/7)

# Equation of the altitude from A
altitude_A = Eq(y, -5)

# Solve for the intersection
intersection = solve((altitude_C, altitude_A), (x, y))
print(intersection)