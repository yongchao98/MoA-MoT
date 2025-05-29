import math

# Coordinates of points A, B, and C
A = (-7, -10)
B = (-2, -3)
C = (-3, -6)

# Vectors AB and BC
AB = (B[0] - A[0], B[1] - A[1])
BC = (C[0] - B[0], C[1] - B[1])

# Dot product of AB and BC
dot_product = AB[0] * BC[0] + AB[1] * BC[1]

# Magnitudes of AB and BC
magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of angle ABC
cos_angle_ABC = dot_product / (magnitude_AB * magnitude_BC)

# Angle ABC in radians
angle_ABC_radians = math.acos(cos_angle_ABC)

# Convert angle to degrees
angle_ABC_degrees = math.degrees(angle_ABC_radians)

# Round the result to 3 decimal places
angle_ABC_degrees_rounded = round(angle_ABC_degrees, 3)

print(angle_ABC_degrees_rounded)