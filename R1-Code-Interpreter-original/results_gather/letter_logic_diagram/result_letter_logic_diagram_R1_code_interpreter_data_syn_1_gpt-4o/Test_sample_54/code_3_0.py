import numpy as np

# Initial grid setup
grid = np.array([
    ['', '', '', '', '', 'd', ''],
    ['g', '', '', 'd', '', 'f', ''],
    ['', '', '', '', 'f', 'e', 'a'],
    ['b', '', '', 'f', 'e', '', ''],
    ['c', '', '', '', 'a', 'g', ''],
    ['d', 'f', '', '', 'g', '', ''],
    ['', 'e', 'a', '', 'b', '', '']
])

# Letters to be used
letters = set('abcdefg')

# Determine the letter for the minor diagonal
minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]

# Find a valid letter for the minor diagonal
minor_diagonal_letter = None
for letter in letters:
    valid = True
    for r, c in minor_diagonal_positions:
        if letter in grid[r, :] or letter in grid[:, c]:
            valid = False
            break
    if valid:
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal
for r, c in minor_diagonal_positions:
    grid[r, c] = minor_diagonal_letter

# Backtracking function to fill the grid
def fill_grid(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return fill_grid(grid, row + 1, 0)
    if grid[row, col] != '':
        return fill_grid(grid, row, col + 1)
    
    for letter in letters:
        if letter not in grid[row, :] and letter not in grid[:, col]:
            grid[row, col] = letter
            if fill_grid(grid, row, col + 1):
                return True
            grid[row, col] = ''
    
    return False

# Fill the grid
fill_grid(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))