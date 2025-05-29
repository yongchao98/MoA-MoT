# Extracting the coordinates from the path data
path_data = "M 22.00,62.00 L 46.00,65.00 L 64.00,60.00 L 91.00,42.00 L 92.00,24.00 L 46.00,19.00 L 22.00,62.00"
# Splitting the path data into commands and coordinates
commands = path_data.split()
# Extracting only the coordinates
coordinates = commands[1::2]
# Converting coordinates to tuples of floats
points = [tuple(map(float, coord.split(','))) for coord in coordinates]
# Using a set to find unique points
unique_points = set(points)
# Counting the number of unique points
num_sides = len(unique_points)
print(num_sides)