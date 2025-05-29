# Initialize the grid
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '0'],
    ['*', '*', '*', '*']
]

# Place the black piece at (3,4)
grid[2][3] = '0'

# Directions: right, left, down, up, diagonals
directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]

# Function to check and flip pieces
def check_and_flip(x, y, dx, dy):
    i, j = x + dx, y + dy
    to_flip = []
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

# Check all directions from the new piece
for dx, dy in directions:
    check_and_flip(2, 3, dx, dy)

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)