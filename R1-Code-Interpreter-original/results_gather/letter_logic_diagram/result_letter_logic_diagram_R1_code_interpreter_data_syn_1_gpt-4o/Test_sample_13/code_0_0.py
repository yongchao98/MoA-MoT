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
# Check which letter can be placed on the minor diagonal
letters = set('abcdefg')
for i in range(7):
    for j in range(7):
        if grid[i, j] in letters:
            letters.remove(grid[i, j])

# Choose a letter for the minor diagonal
minor_diagonal_letter = letters.pop()

# Fill the minor diagonal
for i in range(7):
    grid[i, 6-i] = minor_diagonal_letter

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