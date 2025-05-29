# Define the points from the SVG path
points = [
    (55.57, 80.69),
    (57.38, 65.80),
    (48.90, 57.46),
    (45.58, 47.78),
    (53.25, 36.07),
    (66.29, 48.90),
    (78.69, 61.09),
    (55.57, 80.69)  # Closing the path
]

# Remove the last point as it is the same as the first one (closing the path)
distinct_points = set(points[:-1])

# Count the number of distinct points
num_sides = len(distinct_points)

# Output the number of sides
print(num_sides)