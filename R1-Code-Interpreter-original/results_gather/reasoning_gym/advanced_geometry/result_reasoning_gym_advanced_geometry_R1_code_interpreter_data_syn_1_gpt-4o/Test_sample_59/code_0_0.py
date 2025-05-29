# Define the points
A = (8, -5)
B = (-5, 2)
C = (4, 4)

# Function to calculate slope
def slope(p1, p2):
    if p1[0] == p2[0]:  # vertical line
        return None
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (perpendicular slopes)
slope_altitude_A = -1 / slope_BC if slope_BC is not None else None
slope_altitude_B = -1 / slope_AC if slope_AC is not None else None
slope_altitude_C = -1 / slope_AB if slope_AB is not None else None

# Equation of the altitude from A (using point-slope form)
# y - y1 = m(x - x1)
# y = mx - mx1 + y1
# y = mx + c
if slope_altitude_A is not None:
    c_A = A[1] - slope_altitude_A * A[0]
else:
    c_A = A[0]  # x = constant line

# Equation of the altitude from B
if slope_altitude_B is not None:
    c_B = B[1] - slope_altitude_B * B[0]
else:
    c_B = B[0]  # x = constant line

# Solve the system of equations to find the intersection of altitudes from A and B
# y = slope_altitude_A * x + c_A
# y = slope_altitude_B * x + c_B

if slope_altitude_A is not None and slope_altitude_B is not None:
    # Set the equations equal to each other to find x
    x = (c_B - c_A) / (slope_altitude_A - slope_altitude_B)
    y = slope_altitude_A * x + c_A
elif slope_altitude_A is None:
    x = c_A
    y = slope_altitude_B * x + c_B
elif slope_altitude_B is None:
    x = c_B
    y = slope_altitude_A * x + c_A

# Round the results to 3 decimal places
orthocenter = (round(x, 3), round(y, 3))

print(orthocenter)