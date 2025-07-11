import math

# 1. Define triangle properties
leg_length = 18

# 2. Find the optimal orientation
# We need to maximize k = 2 * floor(Width) + 2 * floor(Height)
# Width = 18 * (cos(phi) + sin(phi))
# Height = 18 * max(cos(phi), sin(phi))
# The maximum is achieved at tan(phi) = 1/2 or tan(phi) = 2.
# Let's use tan(phi) = 1/2
# From tan(phi) = 1/2, we can find sin(phi) and cos(phi)
# sin(phi) = 1 / sqrt(1 + (1/tan(phi))^2) = 1 / sqrt(1 + 2^2) = 1 / sqrt(5)
# cos(phi) = 2 / sqrt(5)
sin_phi = 1 / math.sqrt(5)
cos_phi = 2 / math.sqrt(5)

# 3. Calculate the bounding box dimensions for this orientation
# We assume 0 <= phi <= pi/4, so cos(phi) >= sin(phi)
width = leg_length * (cos_phi + sin_phi)
height = leg_length * cos_phi

# 4. Determine the number of grid lines crossed
# With optimal placement, we can make the number of lines equal to the floor of the dimension
Nv = math.floor(width)
Nh = math.floor(height)

# 5. Calculate the total number of squares
# k = 2 * Nv + 2 * Nh because each line is crossed twice by the convex perimeter
k = 2 * Nv + 2 * Nh

print(f"Optimal orientation is given by tan(phi) = 1/2.")
print(f"Bounding box width = {leg_length} * (2/sqrt(5) + 1/sqrt(5)) = {width:.3f}")
print(f"Bounding box height = {leg_length} * (2/sqrt(5)) = {height:.3f}")
print(f"Number of vertical lines crossed (Nv) = floor({width:.3f}) = {Nv}")
print(f"Number of horizontal lines crossed (Nh) = floor({height:.3f}) = {Nh}")
print(f"Total squares crossed k = 2 * Nv + 2 * Nh = 2 * {Nv} + 2 * {Nh} = {k}")
print(f"\nThe largest number k is {k}.")
