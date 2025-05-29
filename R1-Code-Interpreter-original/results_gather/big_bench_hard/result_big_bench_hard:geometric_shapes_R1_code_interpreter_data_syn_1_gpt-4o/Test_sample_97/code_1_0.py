# Extract points from the path data
path_data = "M 88.00,67.00 L 70.00,18.00 L 38.00,60.00 L 73.00,48.00 L 60.00,62.00 L 88.00,67.00"
# Split the path data by spaces and commas
commands = path_data.replace(',', ' ').split()
points = []

# Iterate over the commands and extract points
for i in range(1, len(commands), 3):
    x, y = float(commands[i]), float(commands[i+1])
    points.append((x, y))

# Remove duplicate points to find unique vertices
unique_points = set(points)

# Output the number of unique points
print(len(unique_points))