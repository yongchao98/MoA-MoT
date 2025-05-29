import math

# Coordinates of the points
A = (3, 3)
B = (-7, -8)
C = (1, 4)

# Calculate the lengths of the sides
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
AC = math.sqrt((C[0] - A[0])**2 + (C[1] - A[1])**2)

# Calculate the angle at B using the law of cosines
cos_B = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)
angle_B_rad = math.acos(cos_B)
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)