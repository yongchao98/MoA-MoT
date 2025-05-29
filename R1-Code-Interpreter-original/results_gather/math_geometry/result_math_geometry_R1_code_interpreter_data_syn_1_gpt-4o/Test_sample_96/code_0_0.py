import math

# Radii of the circles
r1 = 2  # Radius of the first smaller circle
r2 = 3  # Radius of the second smaller circle
R = 5   # Radius of the larger circle

# Calculate the areas
area_large_circle = math.pi * R**2
area_small_circle1 = math.pi * r1**2
area_small_circle2 = math.pi * r2**2

# Calculate the shaded area
shaded_area = area_large_circle - (area_small_circle1 + area_small_circle2)

# Output the shaded area
print(shaded_area)