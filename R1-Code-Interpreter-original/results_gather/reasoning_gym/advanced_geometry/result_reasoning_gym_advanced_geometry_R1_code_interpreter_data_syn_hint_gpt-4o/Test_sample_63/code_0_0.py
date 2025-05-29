import math

# Coordinates
A = (6, 7)
B = (-6, 5)
C = (-10, -6)

# Vectors BA and BC
BA = (A[0] - B[0], A[1] - B[1])
BC = (C[0] - B[0], C[1] - B[1])

# Dot product of BA and BC
dot_product = BA[0] * BC[0] + BA[1] * BC[1]

# Magnitudes of BA and BC
magnitude_BA = math.sqrt(BA[0]**2 + BA[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of angle ABC
cos_angle_ABC = dot_product / (magnitude_BA * magnitude_BC)

# Angle ABC in radians
angle_ABC_radians = math.acos(cos_angle_ABC)

# Convert angle to degrees
angle_ABC_degrees = math.degrees(angle_ABC_radians)

# Round to 3 decimal places
angle_ABC_degrees_rounded = round(angle_ABC_degrees, 3)

print(angle_ABC_degrees_rounded)