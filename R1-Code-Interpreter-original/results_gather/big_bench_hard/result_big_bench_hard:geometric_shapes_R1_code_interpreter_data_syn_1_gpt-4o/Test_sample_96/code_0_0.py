import math

# Define the points
points = [
    (38.00, 62.00),
    (48.00, 60.00),
    (51.00, 49.00),
    (54.00, 60.00),
    (65.00, 62.00),
    (54.00, 64.00),
    (51.00, 74.00),
    (48.00, 64.00),
    (38.00, 62.00)
]

# Function to calculate distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate distances between consecutive points
distances = [distance(points[i], points[i+1]) for i in range(len(points) - 1)]

# Print the distances
print(distances)