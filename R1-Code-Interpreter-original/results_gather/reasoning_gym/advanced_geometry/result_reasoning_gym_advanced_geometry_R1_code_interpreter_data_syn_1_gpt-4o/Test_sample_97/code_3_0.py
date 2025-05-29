import math

# Coordinates of the vertices
A = (6, -2)
B = (-9, 2)
C = (10, -3)

# Function to calculate the distance between two points
def distance(point1, point2):
    return math.sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2)

# Calculate the lengths of the sides
AB = distance(A, B)
BC = distance(B, C)
AC = distance(A, C)

# Use the law of cosines to find the angle at B
# cos(B) = (AB^2 + BC^2 - AC^2) / (2 * AB * BC)
cos_B = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)

# Calculate the angle in radians and then convert to degrees
angle_B_rad = math.acos(cos_B)
angle_B_deg = math.degrees(angle_B_rad)

# Round the angle to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)