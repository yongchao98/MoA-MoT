# Define the points from the path data
points = [
    (49.47, 26.27),
    (55.28, 65.93),
    (48.51, 77.47),
    (34.78, 81.76),
    (36.76, 67.00),
    (14.38, 76.83),
    (49.47, 26.27)  # Closing the path
]

# Use a set to find unique points
unique_points = set(points)

# Print the number of unique points
print(len(unique_points))