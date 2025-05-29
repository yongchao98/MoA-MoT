import math

# Coordinates of the points
A = (-10, 5)
B = (10, 6)
C = (-1, 7)

# Vectors BA and BC
BA = (A[0] - B[0], A[1] - B[1])
BC = (C[0] - B[0], C[1] - B[1])

# Dot product of BA and BC
dot_product = BA[0] * BC[0] + BA[1] * BC[1]

# Magnitudes of BA and BC
magnitude_BA = math.sqrt(BA[0]**2 + BA[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of the angle at B
cos_angle_B = dot_product / (magnitude_BA * magnitude_BC)

# Angle in radians
angle_B_rad = math.acos(cos_angle_B)

# Convert the angle to degrees
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)