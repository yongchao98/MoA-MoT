import math

# Coordinates of the dilated vertices
vertices = [(18, -14), (14, -14), (18, -18), (14, -18)]

# Function to calculate distance from origin
def distance_from_origin(x, y):
    return math.sqrt(x**2 + y**2)

# Calculate distances
distances = [distance_from_origin(x, y) for x, y in vertices]

# Find the index of the maximum distance
max_index = distances.index(max(distances))

# Get the vertex with the maximum distance
farthest_vertex = vertices[max_index]

print(farthest_vertex)