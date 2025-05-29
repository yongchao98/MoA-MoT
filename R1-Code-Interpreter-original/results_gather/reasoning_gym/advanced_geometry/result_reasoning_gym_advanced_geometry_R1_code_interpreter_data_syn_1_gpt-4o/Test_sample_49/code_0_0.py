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

# Calculate the angle in radians
theta_radians = math.acos(cos_theta)

# Convert the angle to degrees
theta_degrees = math.degrees(theta_radians)

# Round the result to 3 decimal places
theta_degrees_rounded = round(theta_degrees, 3)

print(theta_degrees_rounded)