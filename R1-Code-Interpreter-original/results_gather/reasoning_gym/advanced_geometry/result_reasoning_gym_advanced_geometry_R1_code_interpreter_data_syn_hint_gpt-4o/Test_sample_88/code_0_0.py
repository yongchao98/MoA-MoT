import math

# Coordinates of the points
A = (-3, -4)
B = (10, 8)
C = (4, -8)

# Function to calculate distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate side lengths
AB = distance(A, B)
BC = distance(B, C)
AC = distance(A, C)

# Use the law of cosines to find angle at B
# cos(B) = (AB^2 + BC^2 - AC^2) / (2 * AB * BC)
cos_B = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)

# Calculate angle B in degrees
angle_B = math.degrees(math.acos(cos_B))

# Round the angle to 3 decimal places
angle_B_rounded = round(angle_B, 3)

print(angle_B_rounded)