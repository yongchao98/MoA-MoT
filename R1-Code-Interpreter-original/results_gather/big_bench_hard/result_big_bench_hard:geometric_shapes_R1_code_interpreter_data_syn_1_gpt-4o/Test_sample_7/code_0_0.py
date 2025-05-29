# Extracting unique points from the SVG path data
path_data = "M 4.12,72.94 L 38.85,35.38 M 38.85,35.38 L 50.28,77.12 L 30.89,80.59 M 30.89,80.59 L 31.43,90.13 L 22.93,92.01 M 22.93,92.01 L 4.12,72.94"

# Split the path data into commands and coordinates
commands = path_data.split()

# Extract coordinates
coordinates = []
for i in range(1, len(commands), 3):
    x, y = map(float, commands[i].split(','))
    coordinates.append((x, y))

# Get unique points
unique_points = set(coordinates)

# Output the number of unique points
print(len(unique_points))