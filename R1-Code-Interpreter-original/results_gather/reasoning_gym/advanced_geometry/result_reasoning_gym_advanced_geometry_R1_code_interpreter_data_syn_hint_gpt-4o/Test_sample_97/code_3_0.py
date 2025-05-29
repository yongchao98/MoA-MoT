import math

# Vectors
BA = (15, -4)
BC = (19, -5)

# Dot product
dot_product = BA[0] * BC[0] + BA[1] * BC[1]

# Magnitudes
magnitude_BA = math.sqrt(BA[0]**2 + BA[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Cosine of angle B
cos_B = dot_product / (magnitude_BA * magnitude_BC)

# Angle in radians and then in degrees
angle_B_rad = math.acos(cos_B)
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)