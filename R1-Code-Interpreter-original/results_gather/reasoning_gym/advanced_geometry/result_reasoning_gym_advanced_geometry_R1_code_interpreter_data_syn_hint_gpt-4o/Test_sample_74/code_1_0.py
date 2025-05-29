import math

# Vectors AB and BC
AB = (2, 16)
BC = (-2, -15)

# Dot product of AB and BC
dot_product = AB[0] * BC[0] + AB[1] * BC[1]

# Magnitudes of AB and BC
magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of the angle at B
cos_theta = dot_product / (magnitude_AB * magnitude_BC)

# Angle in radians
theta_radians = math.acos(cos_theta)

# Convert angle to degrees
theta_degrees = math.degrees(theta_radians)

# Round to 3 decimal places
theta_degrees_rounded = round(theta_degrees, 3)

print(theta_degrees_rounded)