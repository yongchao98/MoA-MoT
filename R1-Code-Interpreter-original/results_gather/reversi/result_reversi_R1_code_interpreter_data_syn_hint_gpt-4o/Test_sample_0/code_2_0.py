# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Place the black piece at (4,3)
grid[3][2] = '0'

# Function to flip pieces in a given direction
def flip_pieces(grid, start_row, start_col, delta_row, delta_col):
    row, col = start_row + delta_row, start_col + delta_col
    pieces_to_flip = []
    while 0 <= row < 4 and 0 <= col < 4:
        if grid[row][col] == '*':
            break
        if grid[row][col] == '0':
            for r, c in pieces_to_flip:
                grid[r][c] = '0'
            break
        pieces_to_flip.append((row, col))
        row += delta_row
        col += delta_col

# Directions to check: (delta_row, delta_col)
directions = [
    (-1, 0), (1, 0),  # Vertical
    (0, -1), (0, 1),  # Horizontal
    (-1, -1), (-1, 1), (1, -1), (1, 1)  # Diagonal
]

# Apply flipping rules
for delta_row, delta_col in directions:
    flip_pieces(grid, 3, 2, delta_row, delta_col)

# Flatten the grid to the required output format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)