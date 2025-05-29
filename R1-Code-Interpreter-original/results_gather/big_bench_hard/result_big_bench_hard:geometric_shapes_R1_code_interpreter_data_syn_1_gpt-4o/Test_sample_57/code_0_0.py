# Define the points from the SVG path
points = [
    (41.10, 39.10),
    (37.66, 20.53),
    (49.73, 18.29),
    (59.87, 21.21),
    (65.24, 34.62),
    (49.95, 55.71),
    (29.03, 41.34),
    (41.10, 39.10)  # Closing the path
]

# Use a set to find unique points
unique_points = set(points)

# Output the number of unique points
print(len(unique_points))