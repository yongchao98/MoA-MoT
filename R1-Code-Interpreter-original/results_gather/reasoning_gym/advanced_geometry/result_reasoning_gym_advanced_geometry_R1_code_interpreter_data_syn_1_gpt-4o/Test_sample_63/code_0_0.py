import math

# Coordinates
A = (6, 7)
B = (-6, 5)
C = (-10, -6)

# Vectors
AB = (B[0] - A[0], B[1] - A[1])
BC = (C[0] - B[0], C[1] - B[1])

# Dot product
dot_product = AB[0] * BC[0] + AB[1] * BC[1]

# Magnitudes
magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of angle
cos_theta = dot_product / (magnitude_AB * magnitude_BC)

# Angle in radians
theta_radians = math.acos(cos_theta)

# Convert to degrees
theta_degrees = math.degrees(theta_radians)

# Round to 3 decimal places
theta_degrees_rounded = round(theta_degrees, 3)

print(theta_degrees_rounded)