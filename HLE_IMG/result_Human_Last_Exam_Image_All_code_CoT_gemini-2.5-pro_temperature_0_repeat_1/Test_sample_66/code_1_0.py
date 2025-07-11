import math

# Define the given values from the image
radius = 34
tangent_length = 17

# Step 1: Calculate the area of the quadrilateral.
# The quadrilateral is composed of two identical right-angled triangles.
# The legs of each triangle are the radius and the tangent line.
area_of_one_triangle = 0.5 * radius * tangent_length
area_of_quadrilateral = 2 * area_of_one_triangle

# Step 2: Calculate the area of the circular sector.
# First, find the angle of the sector in radians.
# The angle for one triangle is alpha = arctan(opposite/adjacent).
# The total angle is 2 * alpha.
total_angle_rad = 2 * math.atan(tangent_length / radius)

# The area of the sector is 0.5 * r^2 * theta.
area_of_sector = 0.5 * (radius**2) * total_angle_rad

# Step 3: Calculate the total shaded area.
total_shaded_area = area_of_quadrilateral - area_of_sector

# Print the final equation with all the numbers
print("The total shaded area is the area of the quadrilateral minus the area of the circular sector.")
print(f"Final Equation: Area = (2 * (1/2 * {radius} * {tangent_length})) - (1/2 * {radius}^2 * (2 * arctan({tangent_length}/{radius})))")
print(f"Area = {int(area_of_quadrilateral)} - {radius**2} * arctan({tangent_length/radius})")
print(f"Area = {int(area_of_quadrilateral)} - {area_of_sector:.4f}")
print(f"Area = {total_shaded_area:.4f}")