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
            if pieces_to_flip:
                for r, c in pieces_to_flip:
                    grid[r][c] = '0'
            break
        pieces_to_flip.append((row, col))
        row += delta_row
        col += delta_col

# Place the black piece at (1, 2) which is (0, 1) in 0-based index
grid[0][1] = '0'

# Directions: (delta_row, delta_col)
directions = [
    (0, 1),  # right
    (0, -1), # left
    (1, 0),  # down
    (-1, 0), # up
    (1, 1),  # down-right
    (-1, -1),# up-left
    (1, -1), # down-left
    (-1, 1)  # up-right
]

# Check all directions for flipping
for delta_row, delta_col in directions:
    flip_pieces(grid, 0, 1, delta_row, delta_col)

# Flatten the grid to the required format
flattened_grid = [cell for row in grid for cell in row]
print(flattened_grid)