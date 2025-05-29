# Define the points from the SVG path
points = [
    (82.59, 29.92),
    (82.20, 49.38),
    (65.38, 41.46),
    (64.44, 22.89),
    (82.59, 29.92)  # Closing the path
]

# Remove duplicate points to count unique vertices
unique_points = set(points)

# Count the number of unique points
num_sides = len(unique_points)

# Determine the shape based on the number of sides
shape = ""
if num_sides == 3:
    shape = "triangle"
elif num_sides == 4:
    shape = "kite"
elif num_sides == 5:
    shape = "pentagon"
elif num_sides == 6:
    shape = "hexagon"
elif num_sides == 7:
    shape = "heptagon"
elif num_sides == 8:
    shape = "octagon"
else:
    shape = "unknown"

print(shape)