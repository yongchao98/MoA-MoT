import math

# Function to calculate the area of a polygon using the shoelace formula
def polygon_area(vertices):
    n = len(vertices)
    area = 0
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]
        area += x1 * y2 - y1 * x2
    return abs(area) / 2

# Regular octagon side length
s = 1  # Assume side length is 1 for simplicity

# Calculate the radius of the circumscribed circle
R = s / (2 * math.sin(math.pi / 8))

# Calculate coordinates of the vertices of the octagon
vertices = [(R * math.cos(2 * math.pi * i / 8), R * math.sin(2 * math.pi * i / 8)) for i in range(8)]

# Assign vertices to points A, B, C, D, E, F, G, H
A, B, C, D, E, F, G, H = vertices

# Midpoints M and N
M = ((B[0] + C[0]) / 2, (B[1] + C[1]) / 2)
N = ((F[0] + G[0]) / 2, (F[1] + G[1]) / 2)

# Center O
O = (0, 0)

# Calculate areas
area_ABMO = polygon_area([A, B, M, O])
area_EDCMO = polygon_area([E, D, C, M, O])

# Calculate the ratio
ratio = area_ABMO / area_EDCMO

# Output the result
print(f"<<<{ratio}>>>")