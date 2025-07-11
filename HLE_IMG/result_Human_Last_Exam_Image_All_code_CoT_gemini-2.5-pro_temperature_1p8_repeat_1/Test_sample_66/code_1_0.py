import math

# Define the given dimensions
radius = 34
tangent_segment = 17

# Step 1: Calculate the area of the two large triangles formed by the tangents.
# The area of one triangle is (1/2) * base * height = (1/2) * radius * tangent_segment.
# The total area of the two triangles (the quadrilateral) is 2 * (1/2) * radius * tangent_segment.
area_quad = 2 * (0.5 * radius * tangent_segment)

# Step 2: Calculate the area of the total circular sector.
# The half angle of the sector (2*theta in the diagram) can be found from the triangle.
# tan(2*theta) = tangent_segment / radius
# The angle 2*theta in radians is arctan(tangent_segment / radius).
half_angle_rad = math.atan(tangent_segment / radius)

# The total area is the sum of two identical shaded regions.
# Area of one shaded region = Area of one triangle - Area of its corresponding sector
# Area of one triangle = 0.5 * radius * tangent_segment
area_triangle = 0.5 * radius * tangent_segment
# Area of one sector = (1/2) * radius^2 * (2*theta)
# The interpretation from the drawing shows the blue region associated with a 2*theta sector
area_sector_one_side = 0.5 * radius**2 * half_angle_rad
area_one_region = area_triangle - area_sector_one_side

# The total shaded area is the sum of the two identical regions (blue and red).
total_shaded_area = 2 * area_one_region

# Final Equation: Total Area = 2 * (Area of one Triangle - Area of one Sector)
# Total Area = 2 * ( (1/2 * 34 * 17) - (1/2 * 34^2 * arctan(17/34)) )
# Total Area = (2 * 1/2 * 34 * 17) - (2 * 1/2 * 34^2 * arctan(17/34))
# Total Area = (34 * 17) - (34^2 * arctan(0.5))

print("The calculation for the total shaded area is based on the following formula:")
print("Total Area = (2 * (1/2 * base * height)) - (2 * (1/2 * radius^2 * angle_2theta))")
print("Total Area = (radius * height) - (radius^2 * angle_2theta)")
print(f"Total Area = ({radius} * {tangent_segment}) - ({radius}^2 * arctan({tangent_segment}/{radius}))")
print(f"Total Area = {area_quad} - {radius**2} * arctan(0.5)")
print(f"Total Area = {area_quad} - {radius**2 * half_angle_rad:.4f}")
print(f"Total Area = {total_shaded_area:.4f}")
