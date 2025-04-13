import math

# Original vertices
vertices = [(9, -7), (7, -7), (9, -9), (7, -9)]

# Dilation center and scale factor
dilation_center = (0, 0)
scale_factor = 2

# Function to dilate a point
def dilate_point(point, center, scale):
    return (scale * (point[0] - center[0]), scale * (point[1] - center[1]))

# Dilate all vertices
dilated_vertices = [dilate_point(v, dilation_center, scale_factor) for v in vertices]

# Function to calculate distance from origin
def distance_from_origin(point):
    return math.sqrt(point[0]**2 + point[1]**2)

# Find the vertex farthest from the origin
farthest_vertex = max(dilated_vertices, key=distance_from_origin)

# Output the result
print(f"<<<{farthest_vertex}>>>")