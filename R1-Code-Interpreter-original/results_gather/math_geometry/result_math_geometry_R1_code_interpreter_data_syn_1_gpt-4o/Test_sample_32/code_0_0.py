from sympy import symbols, Eq, solve

# Define variables
x, y = symbols('x y')

# Equations of lines AF and BG
eq1 = Eq(y, 3*x)
eq2 = Eq(y, -3/2*x + 15/2)

# Solve for intersection point E
intersection = solve((eq1, eq2), (x, y))
E_x, E_y = intersection[x], intersection[y]

# Coordinates of A, B, and E
A = (0, 0)
B = (5, 0)
E = (E_x, E_y)

# Calculate area of triangle AEB
area = abs(A[0]*(B[1]-E[1]) + B[0]*(E[1]-A[1]) + E[0]*(A[1]-B[1])) / 2
print(area)