import math

# Coordinates of the points
A = (-10, 5)
B = (10, 6)
C = (-1, 7)

# Calculate the lengths of the sides
AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
AC = math.sqrt((C[0] - A[0])**2 + (C[1] - A[1])**2)

# Calculate the cosine of angle ABC using the law of cosines
cos_ABC = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)

# Calculate the angle in degrees
angle_ABC = math.degrees(math.acos(cos_ABC))

# Round the angle to 3 decimal places
angle_ABC_rounded = round(angle_ABC, 3)

print(angle_ABC_rounded)