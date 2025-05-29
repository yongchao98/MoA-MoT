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
# Check which letter can be placed in all diagonal positions
minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
possible_letters = letters.copy()

for r, c in minor_diagonal_positions:
    row_letters = set(grid[r, :])
    col_letters = set(grid[:, c])
    possible_letters &= (letters - row_letters - col_letters)

# Choose the first possible letter for the minor diagonal
minor_diagonal_letter = possible_letters.pop()

# Fill the minor diagonal
for r, c in minor_diagonal_positions:
    grid[r, c] = minor_diagonal_letter

# Fill the rest of the grid
for r in range(7):
    for c in range(7):
        if grid[r, c] == '':
            row_letters = set(grid[r, :])
            col_letters = set(grid[:, c])
            possible_letters = letters - row_letters - col_letters
            grid[r, c] = possible_letters.pop()

# Print the filled grid
for row in grid:
    print(','.join(row))