# Define the coordinates of the two points
x1, y1 = 1, 6
x2, y2 = 17, 4

# Calculate the slope
m = (y2 - y1) / (x2 - x1)

# Calculate the y-intercept of the line
b = y1 - m * x1

# Create a 20x20 grid
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the grid based on the line equation y = mx + b
for x in range(20):
    y = round(m * x + b)
    if 0 <= y < 20:
        grid[y][x] = '✅'

# Print the grid
for row in grid:
    print(''.join(row))