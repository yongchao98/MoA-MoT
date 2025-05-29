import math

# Define the vectors
AB = (10, 1)
BC = (0, -4)

# Calculate the dot product
dot_product = AB[0] * BC[0] + AB[1] * BC[1]

# Calculate the magnitudes
magnitude_AB = math.sqrt(AB[0]**2 + AB[1]**2)
magnitude_BC = math.sqrt(BC[0]**2 + BC[1]**2)

# Calculate the cosine of the angle
cos_theta = dot_product / (magnitude_AB * magnitude_BC)

# Calculate the angle in radians
theta_radians = math.acos(cos_theta)

# Convert the angle to degrees
theta_degrees = math.degrees(theta_radians)

# Round the result to 3 decimal places
theta_degrees_rounded = round(theta_degrees, 3)

print(theta_degrees_rounded)