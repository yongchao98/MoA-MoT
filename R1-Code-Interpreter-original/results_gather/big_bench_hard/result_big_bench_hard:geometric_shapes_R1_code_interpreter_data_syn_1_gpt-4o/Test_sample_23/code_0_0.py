# Extract points from the SVG path data
path_data = "M 59.64,77.72 L 75.18,56.50 M 75.18,56.50 L 6.90,59.13 M 6.90,59.13 L 22.09,77.44 M 22.09,77.44 L 2.73,94.57 M 2.73,94.57 L 91.78,91.66 M 91.78,91.66 L 59.64,77.72"
points = set()

# Split the path data into commands and coordinates
commands = path_data.split(' ')
for i in range(1, len(commands), 3):
    x, y = float(commands[i].strip(',')), float(commands[i+1])
    points.add((x, y))

# Output the number of unique points
print(len(points))