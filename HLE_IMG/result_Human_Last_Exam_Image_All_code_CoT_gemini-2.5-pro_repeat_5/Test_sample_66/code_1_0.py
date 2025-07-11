import math

# Step 1: Define the given parameters from the image.
R = 34  # Radius of the circular sector
d = 17  # Length of the perpendicular line segment

# Step 2: The angle theta can be found using sin(theta) = d / R.
# We don't need the angle itself, just its sine.
sin_theta = d / R

# Step 3: Calculate the area of one elementary triangle formed by two radii and the angle theta.
# The formula for the area of a triangle is 0.5 * a * b * sin(C).
# Here, a = R, b = R, and C = theta.
# Area_triangle = 0.5 * R * R * sin_theta
# This simplifies to 0.5 * R * R * (d / R) = 0.5 * R * d
area_of_one_triangle = 0.5 * R * d

# Step 4: The total shaded region is composed of four such identical regions.
# Based on the principle of quadrature often used in such geometry puzzles,
# the area of each shaded "lune" is equal to the area of its corresponding triangle.
# Therefore, the total shaded area is 4 times the area of one elementary triangle.
number_of_regions = 4
total_area = number_of_regions * area_of_one_triangle

# Print the breakdown of the final calculation
print("The area of one of the four elementary triangles is (1/2) * radius * perpendicular_distance.")
print(f"Area of one triangle = 0.5 * {R} * {d} = {area_of_one_triangle}")
print("The total shaded area is the sum of the areas of the four symmetric regions.")
print(f"Total Shaded Area = {number_of_regions} * {int(area_of_one_triangle)} = {int(total_area)}")
