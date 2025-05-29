# Define the points from the SVG path
points = [
    (9.80, 40.08),
    (12.27, 18.53),
    (29.17, 6.60),
    (43.19, 30.88),
    (27.01, 47.51),
    (9.80, 40.08)  # Closing the shape
]

# Remove duplicate points to find unique vertices
unique_points = set(points)

# Count the number of unique points
num_unique_points = len(unique_points)

# Determine the shape based on the number of unique points
shape = None
if num_unique_points == 3:
    shape = "triangle"
elif num_unique_points == 4:
    shape = "rectangle"  # or possibly a kite
elif num_unique_points == 5:
    shape = "pentagon"
elif num_unique_points == 6:
    shape = "hexagon"
elif num_unique_points == 7:
    shape = "heptagon"
elif num_unique_points == 8:
    shape = "octagon"

print(shape)