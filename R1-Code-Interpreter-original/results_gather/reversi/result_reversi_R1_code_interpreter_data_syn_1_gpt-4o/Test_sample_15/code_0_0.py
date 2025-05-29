# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Move: place a black piece at (2, 1)
grid[2][1] = '0'

# Directions: (row_change, col_change)
directions = [
    (0, 1),  # right
    (0, -1), # left
    (1, 0),  # down
    (-1, 0), # up
    (1, 1),  # bottom-right
    (1, -1), # bottom-left
    (-1, 1), # top-right
    (-1, -1) # top-left
]

# Function to check and flip pieces in a given direction
def flip_pieces(grid, start_row, start_col, direction):
    row, col = start_row + direction[0], start_col + direction[1]
    pieces_to_flip = []

    while 0 <= row < 4 and 0 <= col < 4:
        if grid[row][col] == '*':
            break
        if grid[row][col] == '0':
            for r, c in pieces_to_flip:
                grid[r][c] = '0'
            break
        pieces_to_flip.append((row, col))
        row += direction[0]
        col += direction[1]

# Apply flipping rules in all directions
for direction in directions:
    flip_pieces(grid, 2, 1, direction)

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)