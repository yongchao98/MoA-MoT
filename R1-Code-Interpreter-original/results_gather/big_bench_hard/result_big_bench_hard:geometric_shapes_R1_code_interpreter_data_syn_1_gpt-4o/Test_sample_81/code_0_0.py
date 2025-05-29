# Define the points from the SVG path
points = [
    (32.73, 47.82),
    (41.38, 48.00),
    (45.88, 39.43),
    (46.35, 49.10),
    (55.09, 52.77),
    (45.61, 52.41),
    (41.30, 60.14),
    (40.64, 51.31),
    (32.73, 47.82)  # Closing the path
]

# Count the number of unique points
unique_points = set(points)

# Output the number of unique points
print(len(unique_points))