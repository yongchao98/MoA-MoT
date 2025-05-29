import math

# Define the points
points = [
    (22.06, 65.18),
    (22.58, 64.52),
    (7.52, 52.61),
    (7.00, 53.26),
    (22.06, 65.18)  # Closing the shape
]

# Function to calculate the distance between two points
def distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

# Calculate the lengths of the sides
side_lengths = [distance(points[i], points[i+1]) for i in range(len(points) - 1)]

# Print the side lengths
print(side_lengths)