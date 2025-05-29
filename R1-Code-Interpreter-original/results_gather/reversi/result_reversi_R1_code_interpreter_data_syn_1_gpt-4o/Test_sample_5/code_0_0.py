# Initial grid setup
grid = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Function to check and flip pieces in a given direction
def flip_pieces(grid, start_row, start_col, delta_row, delta_col):
    n = len(grid)
    row, col = start_row + delta_row, start_col + delta_col
    pieces_to_flip = []

    while 0 <= row < n and 0 <= col < n:
        if grid[row][col] == '*':
            break
        if grid[row][col] == '0':
            for r, c in pieces_to_flip:
                grid[r][c] = '0'
            break
        pieces_to_flip.append((row, col))
        row += delta_row
        col += delta_col

# Place the black piece at (1, 2)
grid[1][2] = '0'

# Check all 8 directions
directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
for delta_row, delta_col in directions:
    flip_pieces(grid, 1, 2, delta_row, delta_col)

# Flatten the grid to a single list for output
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)