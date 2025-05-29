# Define the points from the SVG path
points = [
    (19.24, 16.78),
    (35.66, 38.80),
    (35.35, 47.96),
    (28.47, 55.02),
    (24.85, 45.48),
    (14.57, 58.70),
    (19.24, 16.78)  # Closing the path
]

# Count the number of unique points
unique_points = set(points)
print(len(unique_points))