# Coordinates of the vertices
A = (4, 7)
B = (-1, 5)
C = (5, 2)

# Function to calculate the slope of a line given two points
def slope(p1, p2):
    if p2[0] == p1[0]:  # vertical line
        return None
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = -1 / slope_BC if slope_BC is not None else None
slope_altitude_B = -1 / slope_AC if slope_AC is not None else None
slope_altitude_C = -1 / slope_AB if slope_AB is not None else None

# Function to find the equation of a line given a point and a slope
def line_equation(point, slope):
    if slope is None:  # vertical line
        return (1, 0, -point[0])  # x = point[0]
    else:
        # y = mx + c => c = y - mx
        c = point[1] - slope * point[0]
        return (-slope, 1, -c)  # -mx + y = c

# Equations of the altitudes
altitude_A = line_equation(A, slope_altitude_A)
altitude_B = line_equation(B, slope_altitude_B)

# Function to find the intersection of two lines given in the form ax + by + c = 0
def intersection(line1, line2):
    a1, b1, c1 = line1
    a2, b2, c2 = line2
    determinant = a1 * b2 - a2 * b1
    if determinant == 0:
        return None  # Lines are parallel
    x = (b1 * c2 - b2 * c1) / determinant
    y = (a2 * c1 - a1 * c2) / determinant
    return (x, y)

# Find the orthocenter by intersecting two altitudes
orthocenter = intersection(altitude_A, altitude_B)

# Round the coordinates to 3 decimal places
orthocenter = (round(orthocenter[0], 3), round(orthocenter[1], 3))

print(orthocenter)