# Coordinates of the vertices
A = (-8, 10)
B = (4, -10)
C = (-6, -1)

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
slope_altitude_A = -1 / slope_BC if slope_BC is not None else None
slope_altitude_B = -1 / slope_AC if slope_AC is not None else None
slope_altitude_C = -1 / slope_AB if slope_AB is not None else None

# Function to find the intersection of two lines given a point and a slope
def intersection(p1, m1, p2, m2):
    if m1 is None:  # Vertical line through p1
        x = p1[0]
        y = m2 * (x - p2[0]) + p2[1]
    elif m2 is None:  # Vertical line through p2
        x = p2[0]
        y = m1 * (x - p1[0]) + p1[1]
    else:
        x = (m1 * p1[0] - m2 * p2[0] + p2[1] - p1[1]) / (m1 - m2)
        y = m1 * (x - p1[0]) + p1[1]
    return (x, y)

# Find the orthocenter by intersecting two altitudes
# Using altitudes from A and B
orthocenter = intersection(A, slope_altitude_A, B, slope_altitude_B)

# Round the coordinates to 3 decimal places
orthocenter = (round(orthocenter[0], 3), round(orthocenter[1], 3))

print(orthocenter)