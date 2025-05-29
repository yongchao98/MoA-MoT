# Coordinates of the vertices
A = (-4, 8)
B = (2, 10)
C = (-5, 9)

# Function to calculate the slope of a line given two points
def slope(p1, p2):
    if p2[0] - p1[0] == 0:
        return None  # Vertical line
    return (p2[1] - p1[1]) / (p2[0] - p1[0])

# Slopes of the sides
slope_BC = slope(B, C)
slope_AC = slope(A, C)
slope_AB = slope(A, B)

# Slopes of the altitudes (negative reciprocals)
slope_altitude_A = None if slope_BC is None else -1 / slope_BC
slope_altitude_B = None if slope_AC is None else -1 / slope_AC
slope_altitude_C = None if slope_AB is None else -1 / slope_AB

# Function to find the intersection of two lines given in point-slope form
def intersection(slope1, point1, slope2, point2):
    if slope1 is None:  # Vertical line
        x = point1[0]
        y = slope2 * (x - point2[0]) + point2[1]
    elif slope2 is None:  # Vertical line
        x = point2[0]
        y = slope1 * (x - point1[0]) + point1[1]
    else:
        x = (slope1 * point1[0] - slope2 * point2[0] + point2[1] - point1[1]) / (slope1 - slope2)
        y = slope1 * (x - point1[0]) + point1[1]
    return (x, y)

# Find the intersection of the altitudes from A and B
orthocenter = intersection(slope_altitude_A, A, slope_altitude_B, B)

# Round the coordinates to 3 decimal places
orthocenter = (round(orthocenter[0], 3), round(orthocenter[1], 3))

print(orthocenter)