# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Place the black piece at (4,3)
grid[3][2] = '0'

# Directions: right, left, down, up, diagonals
directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]

# Function to check and flip pieces
def flip_pieces(x, y, dx, dy):
    i, j = x + dx, y + dy
    pieces_to_flip = []
    while 0 <= i < 4 and 0 <= j < 4:
        if grid[i][j] == '*':
            break
        if grid[i][j] == '0':
            for px, py in pieces_to_flip:
                grid[px][py] = '0'
            break
        pieces_to_flip.append((i, j))
        i += dx
        j += dy

# Check all directions
for dx, dy in directions:
    flip_pieces(3, 2, dx, dy)

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)