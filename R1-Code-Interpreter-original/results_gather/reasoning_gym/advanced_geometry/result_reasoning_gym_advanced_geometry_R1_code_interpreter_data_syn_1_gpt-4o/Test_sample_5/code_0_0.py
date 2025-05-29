import math

# Define the vectors
u = (0, -18)
v = (-3, -19)

# Calculate the dot product
dot_product = u[0] * v[0] + u[1] * v[1]

# Calculate the magnitudes of the vectors
magnitude_u = math.sqrt(u[0]**2 + u[1]**2)
magnitude_v = math.sqrt(v[0]**2 + v[1]**2)

# Calculate the cosine of the angle
cos_theta = dot_product / (magnitude_u * magnitude_v)

# Calculate the angle in radians
theta_radians = math.acos(cos_theta)

# Convert the angle to degrees
theta_degrees = math.degrees(theta_radians)

# Round the result to 3 decimal places
theta_degrees_rounded = round(theta_degrees, 3)

print(theta_degrees_rounded)