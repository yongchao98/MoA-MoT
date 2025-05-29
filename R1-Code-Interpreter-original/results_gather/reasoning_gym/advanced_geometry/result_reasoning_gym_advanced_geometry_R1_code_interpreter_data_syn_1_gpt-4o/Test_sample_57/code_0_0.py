# Define the points
A = (9, -2)
B = (-6, 1)
C = (-8, -6)

# Function to calculate slope
def slope(p1, p2):
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_CA = slope(C, A)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC
slope_altitude_B = -1 / slope_CA
slope_altitude_C = -1 / slope_AB

# Equations of the altitudes using point-slope form
# Altitude from A
# y - y1 = m(x - x1) => y = m(x - x1) + y1
def altitude_equation(point, slope):
    x1, y1 = point
    return lambda x: slope * (x - x1) + y1

altitude_A = altitude_equation(A, slope_altitude_A)
altitude_B = altitude_equation(B, slope_altitude_B)

# Solve for intersection of altitude_A and altitude_B
# altitude_A(x) = altitude_B(x)
from sympy import symbols, Eq, solve

x = symbols('x')
equation = Eq(altitude_A(x), altitude_B(x))
x_orthocenter = solve(equation)[0]
y_orthocenter = altitude_A(x_orthocenter)

# Round the coordinates to 3 decimal places
x_orthocenter = round(x_orthocenter, 3)
y_orthocenter = round(y_orthocenter, 3)

print((x_orthocenter, y_orthocenter))