# Initial grid setup
grid = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (1, 2) (0-indexed)
grid[1][2] = '0'

# Directions: (dx, dy)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]

# Function to check and flip pieces in a given direction
def flip_in_direction(x, y, dx, dy):
    nx, ny = x + dx, y + dy
    to_flip = []
    while 0 <= nx < 4 and 0 <= ny < 4:
        if grid[nx][ny] == '*':
            break
        if grid[nx][ny] == '0':
            for fx, fy in to_flip:
                grid[fx][fy] = '0'
            break
        to_flip.append((nx, ny))
        nx += dx
        ny += dy

# Check all directions from the placed piece
for dx, dy in directions:
    flip_in_direction(1, 2, dx, dy)

# Flatten the grid to a single list for the final output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)