from sympy import symbols, Eq, solve

# Define the coordinates of the vertices
xA, yA = 2, 6
xB, yB = 1, -6
xC, yC = 5, -10

# Define symbols for the orthocenter coordinates
x, y = symbols('x y')

# Calculate the slopes of the sides
slope_BC = (yC - yB) / (xC - xB)
slope_AC = (yC - yA) / (xC - xA)
slope_AB = (yB - yA) / (xB - xA)

# Calculate the slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_AC
slope_altitude_C = -1 / slope_AB

# Equations of the altitudes
eq1 = Eq(y - yA, slope_altitude_A * (x - xA))
eq2 = Eq(y - yB, slope_altitude_B * (x - xB))

# Solve the system of equations to find the orthocenter
orthocenter = solve((eq1, eq2), (x, y))

# Round the coordinates to 3 decimal places
x_orthocenter = round(orthocenter[x], 3)
y_orthocenter = round(orthocenter[y], 3)

print((x_orthocenter, y_orthocenter))