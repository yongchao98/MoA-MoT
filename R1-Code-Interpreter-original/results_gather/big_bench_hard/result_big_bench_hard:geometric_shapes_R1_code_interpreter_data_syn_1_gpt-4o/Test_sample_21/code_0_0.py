# Define the points from the SVG path
points = [
    (50.55, 76.16),
    (40.00, 72.83),
    (34.22, 59.10),
    (38.51, 46.47),
    (49.40, 60.03),
    (68.54, 44.66),
    (72.08, 58.87),
    (61.72, 73.75),
    (50.55, 76.16)  # Closing the loop
]

# Count unique points
unique_points = set(points)

# Check if the shape is closed and has 6 unique points
is_closed = points[0] == points[-1]
num_unique_points = len(unique_points)

print(is_closed, num_unique_points)