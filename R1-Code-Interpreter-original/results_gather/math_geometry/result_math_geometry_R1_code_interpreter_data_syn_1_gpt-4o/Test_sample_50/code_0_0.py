import math

# Side length of the hexagon
s = 10

# Area of the hexagon
area_hexagon = 6 * (math.sqrt(3) / 4) * s**2

# Area of one equilateral triangle
area_triangle = (math.sqrt(3) / 4) * s**2

# Area of two triangles
area_two_triangles = 2 * area_triangle

# Area of the shaded region
area_shaded = area_hexagon - area_two_triangles

# Output the area of the shaded region in simplest radical form
print(f"{area_shaded / math.sqrt(3)}√3")