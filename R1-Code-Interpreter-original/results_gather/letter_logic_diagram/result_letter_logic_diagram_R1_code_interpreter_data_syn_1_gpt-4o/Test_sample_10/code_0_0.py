import numpy as np

# Initial grid setup
grid = [
    ['d', '', '', 'e', 'b', '', 'g'],
    ['', '', 'e', 'b', '', '', 'd'],
    ['', '', 'b', '', '', '', 'c'],
    ['', '', 'a', '', 'd', 'c', 'f'],
    ['', 'a', 'g', 'd', '', '', 'e'],
    ['', '', '', '', '', '', ''],
    ['', '', '', 'f', '', '', 'a']
]

# Convert grid to numpy array for easier manipulation
grid = np.array(grid)

# Possible letters
letters = set('abcdefg')

# Determine the letter for the minor diagonal
# Check which letter can fit into all diagonal positions
for letter in letters:
    if all(grid[i, 6-i] == '' or grid[i, 6-i] == letter for i in range(7)):
        minor_diagonal_letter = letter
        break

# Fill the minor diagonal
for i in range(7):
    grid[i, 6-i] = minor_diagonal_letter

# Function to fill the grid
def fill_grid(grid):
    for i in range(7):
        for j in range(7):
            if grid[i, j] == '':
                # Determine possible letters for this cell
                row_letters = set(grid[i, :])
                col_letters = set(grid[:, j])
                possible_letters = letters - row_letters - col_letters
                # Assign the first possible letter
                grid[i, j] = possible_letters.pop()
    return grid

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))