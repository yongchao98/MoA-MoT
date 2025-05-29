import math

# Function to calculate the distance between two points
def distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

# Points of the quadrilateral
points = [
    (8.62, 89.32),
    (11.17, 45.94),
    (90.49, 50.60),
    (87.94, 93.98)
]

# Calculate the distances between consecutive points
distances = [
    distance(*points[0], *points[1]),
    distance(*points[1], *points[2]),
    distance(*points[2], *points[3]),
    distance(*points[3], *points[0])
]

# Print the distances
print(distances)