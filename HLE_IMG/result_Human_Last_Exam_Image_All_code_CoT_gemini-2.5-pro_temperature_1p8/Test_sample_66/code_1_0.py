import math

# Define the dimensions from the image
radius = 34
tangent_length = 17

# 1. Calculate the area of one of the large right-angled triangles
area_triangle = 0.5 * radius * tangent_length

# 2. Calculate the angle of the sector in radians
# The angle alpha is given by arctan(opposite/adjacent)
angle_rad = math.atan(tangent_length / radius)

# 3. Calculate the area of the circular sector to be subtracted
area_sector = 0.5 * (radius ** 2) * angle_rad

# 4. Calculate the area of one shaded region
area_one_region = area_triangle - area_sector

# 5. Calculate the total area of the two shaded regions
total_area = 2 * area_one_region

# Print out the explanation and the final equation with all numbers
print("The total shaded area is calculated as 2 * (Area of Triangle - Area of Sector).")
print(f"Area of one triangle = 0.5 * {radius} * {tangent_length} = {area_triangle}")
print(f"Angle of sector in radians = arctan({tangent_length}/{radius}) = {angle_rad:.4f}")
print(f"Area of one sector = 0.5 * {radius}^2 * {angle_rad:.4f} = {area_sector:.4f}")
print(f"Area of one shaded region = {area_triangle} - {area_sector:.4f} = {area_one_region:.4f}")
print(f"Total shaded area = 2 * {area_one_region:.4f} = {total_area:.4f}")

# Final Answer as a full expression
# We can also express the final calculation symbolically:
# Total Area = 2 * (289 - 578 * arctan(0.5))
# Total Area = 578 - 1156 * arctan(0.5)
# Print the final numerical answer
print("\nFinal Answer:")
print(f"Total Area = 2 * ( (1/2) * {radius} * {tangent_length} - (1/2) * {radius}**2 * arctan({tangent_length}/{radius}) )")
print(f"Total Area = {total_area}")