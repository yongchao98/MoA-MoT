import math

# Coordinates of the points
A = (1, 8)
B = (-7, -10)
C = (-1, 6)

# Function to calculate the distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate the lengths of the sides
AB = distance(A, B)
BC = distance(B, C)
AC = distance(A, C)

# Use the law of cosines to find angle at B
# cos(B) = (AB^2 + BC^2 - AC^2) / (2 * AB * BC)
cos_B = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)
angle_B_rad = math.acos(cos_B)

# Convert the angle from radians to degrees
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)