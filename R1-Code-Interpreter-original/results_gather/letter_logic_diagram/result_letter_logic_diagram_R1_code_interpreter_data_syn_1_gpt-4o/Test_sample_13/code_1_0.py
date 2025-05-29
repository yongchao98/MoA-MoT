import numpy as np

# Initial grid setup
grid = [
    ['a', '', 'd', '', 'e', 'c', ''],
    ['b', 'd', 'f', '', '', '', ''],
    ['d', '', '', 'c', '', '', 'a'],
    ['', '', 'c', 'g', '', '', ''],
    ['', '', 'g', '', 'b', '', ''],
    ['', 'g', '', 'b', '', '', 'f'],
    ['', '', '', 'd', '', 'e', '']
]

# Convert grid to numpy array for easier manipulation
grid = np.array(grid)

# Determine the letter for the minor diagonal
# Check if any minor diagonal cells are already filled
minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
minor_diagonal_letter = None

# Check existing letters on the minor diagonal
for i, j in minor_diagonal_indices:
    if grid[i, j] != '':
        minor_diagonal_letter = grid[i, j]
        break

# If no letter is pre-filled on the minor diagonal, choose one
if not minor_diagonal_letter:
    used_letters = set(grid.flatten()) - {''}
    available_letters = set('abcdefg') - used_letters
    minor_diagonal_letter = available_letters.pop()

# Fill the minor diagonal
for i, j in minor_diagonal_indices:
    grid[i, j] = minor_diagonal_letter

# Function to fill the grid
def fill_grid(grid):
    letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if grid[i, j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i, :])
                col_letters = set(grid[:, j])
                possible_letters = letters - row_letters - col_letters
                # Fill the cell with the first possible letter
                grid[i, j] = possible_letters.pop()
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))