import math

# Define the vectors
AB = (2, -12)
BC = (8, 12)

# Calculate the dot product
dot_product = AB[0] * BC[0] + AB[1] * BC[1]

# Calculate the magnitudes
magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Calculate the cosine of the angle
cos_theta = dot_product / (magnitude_AB * magnitude_BC)

# Calculate the angle in radians
angle_radians = math.acos(cos_theta)

# Convert the angle to degrees
angle_degrees = math.degrees(angle_radians)

# Round the angle to 3 decimal places
angle_degrees_rounded = round(angle_degrees, 3)

print(angle_degrees_rounded)