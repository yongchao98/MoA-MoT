import math

# Coordinates of the vertices
A = (4, -4)
B = (-8, 9)
C = (-7, 10)

# Calculate the lengths of the sides
AC = math.sqrt((C[0] - A[0])**2 + (C[1] - A[1])**2)
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)

# Calculate the cosine of angle B using the law of cosines
cos_B = (AC**2 + AB**2 - BC**2) / (2 * AC * AB)

# Calculate the angle in degrees
angle_B = math.degrees(math.acos(cos_B))

# Round the result to 3 decimal places
angle_B_rounded = round(angle_B, 3)

print(angle_B_rounded)