# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '0'],
    ['*', '*', '*', '*']
]

# Place the black piece at (3,4)
grid[2][3] = '0'

# Function to flip pieces in a given direction
def flip_pieces(grid, start, direction, player):
    x, y = start
    dx, dy = direction
    to_flip = []
    x += dx
    y += dy
    while 0 <= x < 4 and 0 <= y < 4:
        if grid[x][y] == '*':
            break
        if grid[x][y] == player:
            for fx, fy in to_flip:
                grid[fx][fy] = player
            break
        to_flip.append((x, y))
        x += dx
        y += dy

# Directions: (dx, dy)
directions = [
    (0, 1),  # Right
    (0, -1), # Left
    (1, 0),  # Down
    (-1, 0), # Up
    (1, 1),  # Down-Right
    (-1, -1),# Up-Left
    (1, -1), # Down-Left
    (-1, 1)  # Up-Right
]

# Check all directions for flipping
for direction in directions:
    flip_pieces(grid, (2, 3), direction, '0')

# Flatten the grid to a single list
flattened_grid = [cell for row in grid for cell in row]

# Print the result
print(flattened_grid)