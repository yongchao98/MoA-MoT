# Initial grid configuration
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Place the black piece at (4, 3)
grid[3][2] = '0'

# Directions: (row_change, col_change)
directions = [
    (-1, 0), (1, 0), (0, -1), (0, 1),  # Up, Down, Left, Right
    (-1, -1), (-1, 1), (1, -1), (1, 1) # Diagonals
]

# Function to check and flip pieces in a given direction
def flip_pieces(grid, start_row, start_col, direction):
    row, col = start_row + direction[0], start_col + direction[1]
    pieces_to_flip = []

    while 0 <= row < 4 and 0 <= col < 4:
        if grid[row][col] == '1':
            pieces_to_flip.append((row, col))
        elif grid[row][col] == '0':
            for r, c in pieces_to_flip:
                grid[r][c] = '0'
            break
        else:
            break
        row += direction[0]
        col += direction[1]

# Check all directions for possible flips
for direction in directions:
    flip_pieces(grid, 3, 2, direction)

# Flatten the grid to the required output format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)