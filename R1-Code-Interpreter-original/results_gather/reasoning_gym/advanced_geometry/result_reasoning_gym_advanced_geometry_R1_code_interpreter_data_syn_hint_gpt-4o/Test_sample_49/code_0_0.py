import math

# Define the vectors
BA = (-3, -7)
BC = (-5, -2)

# Calculate the dot product
dot_product = BA[0] * BC[0] + BA[1] * BC[1]

# Calculate the magnitudes
magnitude_BA = math.sqrt(BA[0]**2 + BA[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Calculate the cosine of the angle
cos_theta = dot_product / (magnitude_BA * magnitude_BC)

# Calculate the angle in degrees
angle_ABC = math.degrees(math.acos(cos_theta))

# Round the angle to 3 decimal places
angle_ABC_rounded = round(angle_ABC, 3)

print(angle_ABC_rounded)