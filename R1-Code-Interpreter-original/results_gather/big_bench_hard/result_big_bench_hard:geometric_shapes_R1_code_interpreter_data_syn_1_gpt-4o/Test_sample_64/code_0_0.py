# Define the points from the SVG path
points = [(93.00, 66.00), (40.00, 65.00), (20.00, 19.00), (93.00, 66.00)]

# Calculate the number of unique points
unique_points = set(points)

# Check if the shape is closed and has three unique points
is_triangle = len(unique_points) == 3 and points[0] == points[-1]

print(is_triangle)