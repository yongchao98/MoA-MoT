# Define the points from the SVG path
points = [
    (10.59, 70.87),
    (29.76, 26.19),
    (73.72, 11.48),
    (88.18, 72.63),
    (42.74, 97.90),
    (10.59, 70.87)  # Closing the path
]

# Count unique points
unique_points = set(points)

# Check if the shape is closed and has 5 unique points
is_closed = points[0] == points[-1]
num_unique_points = len(unique_points)

print(is_closed, num_unique_points)