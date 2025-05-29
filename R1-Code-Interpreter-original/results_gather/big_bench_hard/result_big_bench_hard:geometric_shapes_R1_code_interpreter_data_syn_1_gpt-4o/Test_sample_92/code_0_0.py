# Define the points from the SVG path
points = [(43.79, 3.46), (66.68, 10.86), (56.02, 58.65), (43.79, 3.46)]

# Calculate the number of unique points
unique_points = set(points)

# Check if the shape is closed and has three unique points
is_triangle = len(unique_points) == 3 and points[0] == points[-1]

print(is_triangle)