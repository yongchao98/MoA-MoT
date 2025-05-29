import math

# Coordinates
A = (0, -6)
B = (7, 7)
C = (-1, 6)

# Vectors
AB = (B[0] - A[0], B[1] - A[1])
BC = (C[0] - B[0], C[1] - B[1])

# Dot product
dot_product = AB[0] * BC[0] + AB[1] * BC[1]

# Magnitudes
magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of the angle
cos_theta = dot_product / (magnitude_AB * magnitude_BC)

# Angle in degrees
angle_ABC = math.degrees(math.acos(cos_theta))

# Round to 3 decimal places
angle_ABC = round(angle_ABC, 3)

print(angle_ABC)