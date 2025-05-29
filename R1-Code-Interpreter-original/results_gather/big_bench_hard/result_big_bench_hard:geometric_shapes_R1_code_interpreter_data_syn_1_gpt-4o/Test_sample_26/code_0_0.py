# Define the points from the SVG path
points = [(13.00, 17.00), (78.00, 22.00), (44.00, 59.00), (13.00, 17.00)]

# Find unique points to determine the number of vertices
unique_points = set(points)

# Output the number of unique points
print(len(unique_points))