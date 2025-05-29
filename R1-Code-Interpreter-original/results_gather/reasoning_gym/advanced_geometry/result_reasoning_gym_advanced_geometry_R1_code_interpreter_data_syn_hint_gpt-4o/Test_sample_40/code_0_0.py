from sympy import symbols, Eq, solve

x, y = symbols('x y')
# Equation of the altitude from C
altitude_C = Eq(y - 1, (2/3)*(x - 9))
# Horizontal line y = -5
horizontal_line = Eq(y, -5)

# Solve for intersection
intersection = solve((altitude_C, horizontal_line), (x, y))
print(intersection)