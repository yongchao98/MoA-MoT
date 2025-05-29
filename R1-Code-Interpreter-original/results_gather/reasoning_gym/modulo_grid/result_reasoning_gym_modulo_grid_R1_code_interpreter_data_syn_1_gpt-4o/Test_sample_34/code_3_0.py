# Define the coordinates of the two points
x1, y1 = 1, 6
x2, y2 = 17, 4

# Create a 20x20 grid
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Bresenham's line algorithm
def draw_line(x1, y1, x2, y2):
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx - dy

    while True:
        grid[y1][x1] = '✅'
        if x1 == x2 and y1 == y2:
            break
        e2 = err * 2
        if e2 > -dy:
            err -= dy
            x1 += sx
        if e2 < dx:
            err += dx
            y1 += sy

# Draw the line
draw_line(x1, y1, x2, y2)

# Print the grid
for row in grid:
    print(''.join(row))