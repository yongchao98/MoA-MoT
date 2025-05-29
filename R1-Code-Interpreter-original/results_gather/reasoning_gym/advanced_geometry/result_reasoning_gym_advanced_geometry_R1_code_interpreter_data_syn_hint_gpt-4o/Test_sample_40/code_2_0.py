from sympy import symbols, Eq, solve

x, y = symbols('x y')
# Equation of the altitude from C
altitude_C = Eq(y - 1, (2/3)*(x - 9))
# Horizontal line y = 7
horizontal_line_y7 = Eq(y, 7)

# Solve for intersection
intersection_y7 = solve((altitude_C, horizontal_line_y7), (x, y))
print(intersection_y7)