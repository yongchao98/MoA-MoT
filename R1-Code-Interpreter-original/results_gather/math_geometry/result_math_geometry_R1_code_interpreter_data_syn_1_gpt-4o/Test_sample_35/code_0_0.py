from sympy import symbols, Eq, solve

# Define the variables
x, y = symbols('x y')

# Define the equations of the lines
line1 = Eq(y, 2*x - 4)
line2 = Eq(y, -3*x + 16)

# Find the intersection of the two lines
intersection = solve((line1, line2), (x, y))

# Find the y-intercepts
y_intercept_line1 = solve(line1.subs(x, 0), y)[0]
y_intercept_line2 = solve(line2.subs(x, 0), y)[0]

# Vertices of the triangle
vertices = [(0, y_intercept_line1), (0, y_intercept_line2), (intersection[x], intersection[y])]

# Calculate the area of the triangle using the formula:
# Area = 0.5 * | x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2) |
x1, y1 = vertices[0]
x2, y2 = vertices[1]
x3, y3 = vertices[2]

area = 0.5 * abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))

print(area)