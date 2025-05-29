# Define the points from the SVG path
points = [(73.73, 87.22), (79.29, 59.06), (53.87, 38.89), (73.73, 87.22)]

# Calculate the number of unique points
unique_points = set(points)

# Check if the shape is closed and has three unique points
is_triangle = len(unique_points) == 3 and points[0] == points[-1]

print(is_triangle)