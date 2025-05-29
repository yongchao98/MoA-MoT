import math

# Vectors BA and BC
BA = (-6, 1)
BC = (1, 10)

# Dot product of BA and BC
dot_product = BA[0] * BC[0] + BA[1] * BC[1]

# Magnitudes of BA and BC
magnitude_BA = math.sqrt(BA[0]**2 + BA[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of the angle
cos_theta = dot_product / (magnitude_BA * magnitude_BC)

# Angle in radians
theta_radians = math.acos(cos_theta)

# Convert to degrees
theta_degrees = math.degrees(theta_radians)

# Round to 3 decimal places
theta_degrees_rounded = round(theta_degrees, 3)

print(theta_degrees_rounded)