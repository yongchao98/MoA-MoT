import math

# Side length of the equilateral triangle
s = 4

# Area of the equilateral triangle
area_triangle = (math.sqrt(3) / 4) * s**2

# Radius of the circle
r = 4

# Central angle in radians
theta = math.pi / 3

# Area of one sector
area_sector = 0.5 * r**2 * theta

# Total area of three sectors
total_area_sectors = 3 * area_sector

# Area of the region
area_region = area_triangle - total_area_sectors

# Expressing in the form a√b + cπ
a = 4
b = 3
c = -8

# Calculating a + b + c
result = a + b + c

print(result)