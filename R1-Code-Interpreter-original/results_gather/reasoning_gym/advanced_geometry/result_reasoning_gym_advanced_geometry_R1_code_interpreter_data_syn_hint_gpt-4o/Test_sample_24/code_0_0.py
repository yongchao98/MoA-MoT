import math

# Coordinates of the points
A = (-10, 5)
B = (10, 6)
C = (-1, 7)

# Calculate the lengths of the sides
AC = math.sqrt((A[0] - C[0])**2 + (A[1] - C[1])**2)
AB = math.sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2)
BC = math.sqrt((B[0] - C[0])**2 + (B[1] - C[1])**2)

# Calculate the angle at B using the law of cosines
cos_B = (AC**2 + BC**2 - AB**2) / (2 * AC * BC)
angle_B_rad = math.acos(cos_B)

# Convert the angle from radians to degrees
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)