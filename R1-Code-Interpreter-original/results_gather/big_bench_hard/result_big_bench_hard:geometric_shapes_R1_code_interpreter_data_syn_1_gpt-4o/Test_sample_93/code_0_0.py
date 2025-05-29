import math

# Define the points
points = [
    (36.48, 46.06),
    (35.35, 54.14),
    (73.45, 59.46),
    (74.58, 51.38),
    (36.48, 46.06)  # Closing the shape
]

# Function to calculate the distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate the distances between consecutive points
distances = [distance(points[i], points[i+1]) for i in range(len(points) - 1)]

# Print the distances
print(distances)