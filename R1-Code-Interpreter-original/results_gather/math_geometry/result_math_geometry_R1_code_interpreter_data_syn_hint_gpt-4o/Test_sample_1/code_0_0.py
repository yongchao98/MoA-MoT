import math

# Side length of the octagon
s = 1  # Assume side length is 1 for simplicity

# Area of the regular octagon
area_octagon = 2 * (1 + math.sqrt(2)) * s**2

# Area of one of the triangles formed by the center and a side
area_triangle = 0.5 * s * s * math.sin(math.pi / 4)

# Area of the quadrilateral ABMO
area_ABMO = 2 * area_triangle

# Area of the quadrilateral EDCMO
area_EDCMO = 3 * area_triangle

# Ratio of the areas
ratio = area_ABMO / area_EDCMO

print(ratio)