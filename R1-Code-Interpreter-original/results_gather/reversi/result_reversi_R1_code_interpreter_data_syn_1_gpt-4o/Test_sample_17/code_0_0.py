# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Place the black piece at (4, 3)
grid[3][2] = '0'

# Directions: horizontal, vertical, and diagonal
directions = [
    (0, 1), (0, -1),  # Horizontal
    (1, 0), (-1, 0),  # Vertical
    (1, 1), (-1, -1), # Diagonal \
    (1, -1), (-1, 1)  # Diagonal /
]

# Function to check and flip pieces
def flip_pieces(x, y, dx, dy):
    to_flip = []
    i, j = x + dx, y + dy
    while 0 <= i < 4 and 0 <= j < 4:
        if grid[i][j] == '*':
            break
        if grid[i][j] == '0':
            for fx, fy in to_flip:
                grid[fx][fy] = '0'
            break
        to_flip.append((i, j))
        i += dx
        j += dy

# Check all directions from the placed piece
for dx, dy in directions:
    flip_pieces(3, 2, dx, dy)

# Flatten the grid to the required output format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)