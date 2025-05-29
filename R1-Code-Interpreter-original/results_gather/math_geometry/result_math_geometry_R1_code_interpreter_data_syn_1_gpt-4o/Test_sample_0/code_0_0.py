import math

# Original vertices after dilation
vertices = [(14, -14), (18, -14), (14, -18), (18, -18)]

# Calculate distances from the origin
distances = [(x, y, math.sqrt(x**2 + y**2)) for x, y in vertices]

# Find the vertex with the maximum distance
farthest_vertex = max(distances, key=lambda item: item[2])

# Output the farthest vertex
print(farthest_vertex[:2])