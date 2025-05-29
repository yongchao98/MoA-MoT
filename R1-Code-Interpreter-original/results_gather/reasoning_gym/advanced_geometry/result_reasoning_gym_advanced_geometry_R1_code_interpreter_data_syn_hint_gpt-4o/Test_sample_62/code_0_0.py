import math

# Coordinates of the points
A = (-5, -7)
B = (-2, -1)
C = (3, -10)

# Calculate the lengths of the sides
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
AC = math.sqrt((C[0] - A[0])**2 + (C[1] - A[1])**2)

# Calculate the cosine of the angle at B using the law of cosines
cos_ABC = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)

# Calculate the angle in radians and then convert to degrees
angle_ABC_rad = math.acos(cos_ABC)
angle_ABC_deg = math.degrees(angle_ABC_rad)

# Round the result to 3 decimal places
angle_ABC_deg_rounded = round(angle_ABC_deg, 3)

print(angle_ABC_deg_rounded)