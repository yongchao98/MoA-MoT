import math

# Define the given values
radius = 34
height = 17

# 1. Calculate the area of one right-angled triangle
area_one_triangle = 0.5 * radius * height

# 2. Calculate the total area of the four identical triangles
total_area_triangles = 4 * area_one_triangle

# 3. Calculate the angle theta in radians for one triangle
# tan(theta) = height / radius
theta_rad = math.atan(height / radius)

# 4. Calculate the total angle for the large sector
total_angle_rad = 4 * theta_rad

# 5. Calculate the area of the large circular sector
# Formula: (angle_rad / 2) * radius^2
area_sector = (total_angle_rad / 2) * (radius**2)

# 6. Calculate the total shaded area
total_shaded_area = total_area_triangles - area_sector

# Print the final equation with all the numbers
print("The total shaded area is the area of the four outer triangles minus the area of the inner circular sector.")
print(f"Final Equation: 4 * (0.5 * {radius} * {height}) - 2 * {radius}**2 * arctan({height}/{radius})")
print(f"Calculation: {total_area_triangles} - {area_sector:.4f} = {total_shaded_area:.4f}")
print(f"The total area of the shaded regions is {total_shaded_area:.4f}")